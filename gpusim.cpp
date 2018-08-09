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

#include <sys/time.h>

#include "fingerprintdb_cuda.h"
#include "local_qinfo.h"

using std::shared_ptr;
using std::vector;
using gpusim::FingerprintDB;
using gpusim::Fingerprint;

namespace gpusim
{
GPUSimServer::GPUSimServer(const QStringList& database_fnames)
{
    for(auto database_fname : database_fnames) {
        // Read from .fsim file into byte arrays
        QSize fingerprint_size;
        QByteArray fingerprint_data;
        vector<char*> smiles_vector;
        vector<char*> ids_vector;
        extractData(database_fname, fingerprint_size, fingerprint_data,
                smiles_vector, ids_vector);

        // Create new FingerprintDB for querying on GPU
        auto fps = std::shared_ptr<FingerprintDB>(new FingerprintDB(
                fingerprint_size.width(), fingerprint_size.height(),
                fingerprint_data.data(), smiles_vector, ids_vector));

        QFileInfo file_info(database_fname);
        QString socket_name = file_info.baseName();
        if (!setupSocket(socket_name))
            return;
        qDebug() << "Setting up DB:  " << socket_name;
        m_databases[socket_name] = fps;
    }

    // Now that we know how much total memory is required, divvy it up
    // and allow the fingerprint databases to copy up data
    for(auto db : m_databases) {
        db->copyToGPU(get_available_gpu_memory());
    }

    qInfo() << "Ready for searches.";
};

void GPUSimServer::extractData(const QString& database_fname,
                                QSize& fingerprint_size,
                                QByteArray& fingerprint_data,
                                vector<char*>& smiles_vector,
                                vector<char*>& ids_vector)
{
    QFile file(database_fname);
    file.open(QIODevice::ReadOnly);
    QDataStream datastream(&file);
    // Set version so that files will be usable cross-release
    datastream.setVersion(QDataStream::Qt_5_2);

    QByteArray smi_data, id_data;
    datastream >> fingerprint_size;
    datastream >> fingerprint_data;
    datastream >> smi_data;
    datastream >> id_data;

    int count = fingerprint_size.height();
    smiles_vector.resize(count);
    ids_vector.resize(count);

    // Extract smiles vector from serialized data
    QDataStream smi_stream(smi_data);
    for (int i = 0; i < count; i++) {
        smi_stream >> smiles_vector[i];
    }

    // Extract ID vector from serialized data
    QDataStream id_stream(id_data);
    for (int i = 0; i < count; i++) {
        id_stream >> ids_vector[i];
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
        vector<float>& results_scores, unsigned int return_count, CalcType calc_type)
{
    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, nullptr);

    if(calc_type == CalcType::GPU) {
        m_databases[dbname]->search(reference, results_smiles, results_ids,
                results_scores, return_count);
    } else {
        m_databases[dbname]->search_cpu(reference, results_smiles, results_ids,
                results_scores, return_count);
    }

    gettimeofday(&tval_after, nullptr);

    timersub(&tval_after, &tval_before, &tval_result);

    // Using printf here because decimal formatting like this is a PITA w/ cout
    printf("Search completed, time elapsed: %ld.%06ld\n",
           (long int) tval_result.tv_sec, (long int) tval_result.tv_usec);
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

void GPUSimServer::incomingSearchRequest()
{
    auto clientConnection = static_cast<QLocalSocket*>(sender());
    QString dbname = clientConnection->property("dbname").toString();

    // Read incoming Fingerprint binary and put it in Fingerprint object
    QByteArray data = clientConnection->readAll();
    QDataStream qds(&data, QIODevice::ReadOnly);
    int results_requested;
    qds >> results_requested;
    QByteArray fp_data;
    qds >> fp_data;
    const int* raw_fp_data = reinterpret_cast<const int*>(fp_data.constData());
    const int fp_int_size = fp_data.size() / sizeof(int);

    Fingerprint query(fp_int_size);
    query.assign(raw_fp_data, raw_fp_data + fp_int_size);

    // Perform similarity search and return results in relevant vectors
    vector<char *> results_smiles, results_ids;
    vector<float> results_scores;
    
    similaritySearch(query, dbname, results_smiles, results_ids,
            results_scores, results_requested,
            usingGPU() ? CalcType::GPU : CalcType::CPU);

    // Create QByteArrays and QDataStreams to write to corresponding arrays
    QByteArray output_smiles, output_ids, output_scores;
    QDataStream smiles_stream(&output_smiles, QIODevice::WriteOnly);
    QDataStream ids_stream(&output_ids, QIODevice::WriteOnly);
    QDataStream scores_stream(&output_scores, QIODevice::WriteOnly);
    for (int i = 0; i < results_requested; i++) {
        smiles_stream << results_smiles[i];
        ids_stream << results_ids[i];
        scores_stream << results_scores[i];
    }

    // Transmit binary data to client and flush the buffered data
    QByteArray rcount;
    QDataStream rcount_qds(&rcount, QIODevice::WriteOnly);
    rcount_qds << results_requested;
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
