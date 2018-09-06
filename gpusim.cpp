/* -------------------------------------------------------------------------
 * Implements gpusim::GPUSimServer
 * Reads a binary .fsim file and creates a FingerprintDB on GPU.  Uses a local
 * socket to communicate with querying processes.  Fingerprint agnostic.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "gpusim.h"

#include <QByteArray>
#include <QCoreApplication>
#include <QDataStream>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QLocalServer>
#include <QLocalSocket>
#include <QSize>
#include <QTime>

#include <algorithm>
#include <exception>
#include <math.h>

#include "fingerprintdb_cuda.h"
#include "local_qinfo.h"

using std::pair;
using std::shared_ptr;
using std::vector;
using gpusim::FingerprintDB;
using gpusim::Fingerprint;

const int DATABASE_VERSION = 2;

namespace gpusim
{
GPUSimServer::GPUSimServer(const QStringList& database_fnames, int gpu_bitcount)
{
    qDebug() << "--------------------------";
    qDebug() << "Starting up GPUSim Server";
    qDebug() << "--------------------------";
    for(auto database_fname : database_fnames) {
        // Read from .fsim file into byte arrays
        int fp_bitcount, fp_count;
        vector<vector<char> > fingerprint_data;
        vector<char*> smiles_vector;
        vector<char*> ids_vector;
        qDebug() << "Extracting data";
        extractData(database_fname, fp_bitcount, fp_count, fingerprint_data,
                smiles_vector, ids_vector);
        qDebug() << "Finished extracting data";

        qDebug() << "Utilizing " << get_gpu_count() << " GPUs for calculation.";

        // Create new FingerprintDB for querying on GPU
        auto fps = std::shared_ptr<FingerprintDB>(new FingerprintDB(
                fp_bitcount, fp_count,
                fingerprint_data, smiles_vector, ids_vector));

        QFileInfo file_info(database_fname);
        QString socket_name = file_info.baseName();
        if (!setupSocket(socket_name))
            return;
        qDebug() << "Setting up DB:  " << socket_name;
        m_databases[socket_name] = fps;
    }
    if (!setupSocket("all"))
        return;

    // Now that we know how much total memory is required, divvy it up
    // and allow the fingerprint databases to copy up data
    size_t total_db_memory = 0;
    unsigned int max_compounds_in_db = 0;
    int max_fp_bitcount = 0;
    for(auto db : m_databases) {
        total_db_memory += db->getFingerprintDataSize();
        max_compounds_in_db = std::max(max_compounds_in_db, db->count());
        max_fp_bitcount = std::max(max_fp_bitcount,
                db->getFingerprintBitcount());
    }


    unsigned int fold_factor = 1;
    auto gpu_memory = get_available_gpu_memory();

    // Reserve space for the indices vector during search
    gpu_memory -= sizeof(int) * max_compounds_in_db;

    qDebug() << "Database:  " << total_db_memory/1024/1024 << "MB GPU Memory: "  <<
        gpu_memory/1024/1024 << "MB";

    if(total_db_memory > gpu_memory) 
    {
        fold_factor = ceilf(
                (float)total_db_memory / (float)gpu_memory);
        qDebug() << "Folding databases by at least " << fold_factor << " to fit in gpu memory";
    }

    if(gpu_bitcount > 0) {
        unsigned int arg_fold_factor = max_fp_bitcount / gpu_bitcount;
        if (arg_fold_factor < fold_factor) {
            throw std::invalid_argument("GPU bitset not sufficiently small to fit on GPU");
        } 
        fold_factor = arg_fold_factor;
    }

    qInfo() << "Putting graphics card data up.";
    if(usingGPU())
    {
        for(auto db : m_databases) {
            db->copyToGPU(fold_factor);
        }
    }
    qInfo() << "Finished putting graphics card data up.";

    qInfo() << "Ready for searches.";
};

bool GPUSimServer::usingGPU()
{
    return m_use_gpu && (get_gpu_count() != 0);
}

void GPUSimServer::extractData(const QString& database_fname,
                                int& fp_bitcount,
                                int& fp_count,
                                vector<vector<char> >& fingerprint_data,
                                vector<char*>& smiles_vector,
                                vector<char*>& ids_vector)
{
    qDebug() << "Begin extracting data";
    QFile file(database_fname);
    file.open(QIODevice::ReadOnly);
    QDataStream datastream(&file);
    // Set version so that files will be usable cross-release
    datastream.setVersion(QDataStream::Qt_5_2);
    int version;
    datastream >> version;
    if(version != DATABASE_VERSION) {
        throw std::runtime_error("Database version incompatible with this GPUSim version");
    }

    datastream >> fp_bitcount;
    datastream >> fp_count;

    smiles_vector.reserve(fp_count);
    ids_vector.reserve(fp_count);

    int fp_qba_count;
    datastream >> fp_qba_count;
    fingerprint_data.resize(fp_qba_count);
    for(auto& fp_vector : fingerprint_data) {
        QByteArray fp_qba;
        datastream >> fp_qba;
        fp_vector.reserve(fp_qba.size());
        fp_vector.insert(fp_vector.begin(), fp_qba.data(), fp_qba.data()+fp_qba.size());
    }

    int smi_qba_count;
    datastream >> smi_qba_count;
    for(int i=0; i<smi_qba_count; i++) {
        QByteArray smi_qba;
        datastream >> smi_qba;
        // Extract smiles vector from serialized data
        QDataStream smi_stream(smi_qba);
        while(!smi_stream.atEnd()) {
            char* smi;
            smi_stream >> smi;
            smiles_vector.push_back(smi);
        }
    }

    int id_qba_count;
    datastream >> id_qba_count;
    for(int i=0; i<id_qba_count; i++) {
        QByteArray id_qba;
        datastream >> id_qba;
        // Extract smiles vector from serialized data
        QDataStream id_stream(id_qba);
        while(!id_stream.atEnd()) {
            char* id;
            id_stream >> id;
            ids_vector.push_back(id);
        }
    }
}

bool GPUSimServer::setupSocket(const QString& socket_name)
{
    auto server = new QLocalServer(this);
    if (!server->listen(socket_name)) {
        QString socket_location = QString("/tmp/%1").arg(socket_name);
        QFile::remove(socket_location);
        if (!server->listen(socket_name)) {
            qDebug() << "Server start failed on " << socket_location;
            auto app = QCoreApplication::instance();
            app->exit(1);
            return false;
        }
    }

    QLocalServer::connect(server, &QLocalServer::newConnection, this,
                          &GPUSimServer::newConnection);

    return true;
}

void GPUSimServer::similaritySearch(const Fingerprint& reference,
        const QString& dbname,
        vector<char*>& results_smiles, vector<char*>& results_ids,
        vector<float>& results_scores, unsigned int return_count,
        float similarity_cutoff, CalcType calc_type)
{
    if(calc_type == CalcType::GPU) {
        m_databases[dbname]->search(reference, results_smiles, results_ids,
                results_scores, return_count, similarity_cutoff);
    } else {
        m_databases[dbname]->search_cpu(reference, results_smiles, results_ids,
                results_scores, return_count, similarity_cutoff);
    }
};

void GPUSimServer::newConnection()
{
    auto server = dynamic_cast<QLocalServer*>(sender());
    QLocalSocket* clientConnection = server->nextPendingConnection();

    // HACK:  For some reason QLocalSocket->serverName doesn't work in socket
    clientConnection->setProperty("dbname", server->serverName());

    QObject::connect(clientConnection, &QLocalSocket::disconnected,
                     clientConnection, &QLocalSocket::deleteLater);
    QObject::connect(clientConnection, &QLocalSocket::readyRead, this,
                     &GPUSimServer::incomingSearchRequest);
}

void GPUSimServer::searchAll(const Fingerprint& query, int results_requested,
        float similarity_cutoff, vector<char *>&  results_smiles,
        vector<char *>& results_ids, vector<float>& results_scores)
{
    typedef pair<char*, char*> ResultData;
    typedef pair<float, ResultData > SortableResult;

    vector<SortableResult> sortable_results;
    for(auto local_dbname : m_databases.keys()) {
        vector<char *> l_results_smiles, l_results_ids;
        vector<float> l_results_scores;
        similaritySearch(query, local_dbname, l_results_smiles, l_results_ids,
                l_results_scores, results_requested, similarity_cutoff,
                usingGPU() ? CalcType::GPU : CalcType::CPU);
        for(unsigned int i=0; i<l_results_smiles.size(); i++) {
            sortable_results.push_back(SortableResult(l_results_scores[i], 
                        ResultData(l_results_smiles[i], l_results_ids[i])));
        }
    }
    std::sort(sortable_results.begin(), sortable_results.end());
    std::reverse(sortable_results.begin(), sortable_results.end());

    for(auto result : sortable_results) {
        results_scores.push_back(result.first);
        results_smiles.push_back(result.second.first);
        results_ids.push_back(result.second.second);
    }
    int result_size = std::min(results_requested, (int)results_scores.size());
    results_scores.resize(result_size);
    results_smiles.resize(result_size);
    results_ids.resize(result_size);
}

void GPUSimServer::incomingSearchRequest()
{
    auto clientConnection = static_cast<QLocalSocket*>(sender());
    QString dbname = clientConnection->property("dbname").toString();

    // Read incoming Fingerprint binary and put it in Fingerprint object
    QByteArray data = clientConnection->readAll();
    QDataStream qds(&data, QIODevice::ReadOnly);
    int results_requested;
    qds >> results_requested;
    float similarity_cutoff;
    qds >> similarity_cutoff;
    QByteArray fp_data;
    qds >> fp_data;
    const int* raw_fp_data = reinterpret_cast<const int*>(fp_data.constData());
    const int fp_int_size = fp_data.size() / sizeof(int);

    Fingerprint query(fp_int_size);
    query.assign(raw_fp_data, raw_fp_data + fp_int_size);

    // Perform similarity search and return results in relevant vectors
    vector<char *> results_smiles, results_ids;
    vector<float> results_scores;
    
    QTime timer;
    timer.start();

    if(dbname == "all") {
        searchAll(query, results_requested, similarity_cutoff, results_smiles,
                results_ids, results_scores);
    } else {
        similaritySearch(query, dbname, results_smiles, results_ids,
                results_scores, results_requested, similarity_cutoff,
                usingGPU() ? CalcType::GPU : CalcType::CPU);
    }

    qDebug() << "Search completed, time elapsed: " << (float)timer.elapsed() / 1000.0f;

    // Create QByteArrays and QDataStreams to write to corresponding arrays
    QByteArray output_smiles, output_ids, output_scores;
    QDataStream smiles_stream(&output_smiles, QIODevice::WriteOnly);
    QDataStream ids_stream(&output_ids, QIODevice::WriteOnly);
    QDataStream scores_stream(&output_scores, QIODevice::WriteOnly);
    for (unsigned int i = 0; i < results_smiles.size(); i++) {
        smiles_stream << results_smiles[i];
        ids_stream << results_ids[i];
        scores_stream << results_scores[i];
    }

    // Transmit binary data to client and flush the buffered data
    QByteArray rcount;
    QDataStream rcount_qds(&rcount, QIODevice::WriteOnly);
    rcount_qds << (int)results_smiles.size();
    clientConnection->write(rcount);
    clientConnection->write(output_smiles);
    clientConnection->write(output_ids);
    clientConnection->write(output_scores);
    clientConnection->flush();
}


Fingerprint GPUSimServer::getFingerprint(const int index, const QString& dbname)
{
    return m_databases[dbname]->getFingerprint(index);
}

} // end gpusim
