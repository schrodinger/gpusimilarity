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
#include <QThread>
#include <QThreadPool>
#include <QTime>

#include <algorithm>
#include <exception>
#include <map>
#include <math.h>
#include <set>
#include <sstream>

#include "fingerprintdb_cuda.h"
#include "local_qinfo.h"

using gpusim::Fingerprint;
using gpusim::FingerprintDB;
using std::map;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;

const int DATABASE_VERSION = 3;

namespace gpusim
{

class DecompressAssignFPRunnable : public QRunnable
{
  public:
    QByteArray compressed_data;
    vector<char>& fp_vector;

    DecompressAssignFPRunnable(vector<char>& fp_in) : fp_vector(fp_in){};
    void run() override
    {
        QByteArray fp_qba = qUncompress(compressed_data);
        compressed_data.clear();
        fp_vector.reserve(fp_qba.size());
        fp_vector.insert(fp_vector.begin(), fp_qba.data(),
                         fp_qba.data() + fp_qba.size());
    }
};

class DecompressAssignStringRunnable : public QRunnable
{
  public:
    QByteArray compressed_data;
    vector<char*>& string_vector;

    DecompressAssignStringRunnable(vector<char*>& string_in)
        : string_vector(string_in){};
    void run() override
    {
        QByteArray string_qba = qUncompress(compressed_data);
        compressed_data.clear();
        QDataStream string_stream(string_qba);
        while (!string_stream.atEnd()) {
            char* smi;
            string_stream >> smi;
            string_vector.push_back(smi);
        }
        string_qba.clear();
    }
};

GPUSimServer::GPUSimServer(const QStringList& database_fnames, int gpu_bitcount)
{
    qDebug() << "--------------------------";
    qDebug() << "Starting up GPUSim Server";
    qDebug() << "--------------------------";
    qDebug() << "Utilizing" << get_gpu_count() << "GPUs for calculation.";

    if (!setupSocket())
        return;

    for (auto database_fname : database_fnames) {
        // Read from .fsim file into byte arrays
        int fp_bitcount, fp_count;
        QString dbkey;
        vector<vector<char>> fingerprint_data;
        vector<char*> smiles_vector;
        vector<char*> ids_vector;
        qDebug() << "Extracting data:" << database_fname;
        extractData(database_fname, fp_bitcount, fp_count, dbkey,
                    fingerprint_data, smiles_vector, ids_vector);
        qDebug() << "Finished extracting data";

        // Create new FingerprintDB for querying on GPU
        auto fps = std::shared_ptr<FingerprintDB>(
            new FingerprintDB(fp_bitcount, fp_count, dbkey, fingerprint_data,
                              smiles_vector, ids_vector));

        QFileInfo file_info(database_fname);
        QString db_name = file_info.baseName();
        m_databases[db_name] = fps;
    }

    // Now that we know how much total memory is required, divvy it up
    // and allow the fingerprint databases to copy up data
    size_t total_db_memory = 0;
    unsigned int max_compounds_in_db = 0;
    int max_fp_bitcount = 0;
    for (auto db : m_databases) {
        total_db_memory += db->getFingerprintDataSize();
        max_compounds_in_db = std::max(max_compounds_in_db, db->count());
        max_fp_bitcount =
            std::max(max_fp_bitcount, db->getFingerprintBitcount());
    }

    unsigned int fold_factor = 1;
    auto gpu_memory = get_available_gpu_memory();

    // Reserve space for the indices vector during search
    gpu_memory -= sizeof(int) * max_compounds_in_db;

    qDebug() << "Database:  " << total_db_memory / 1024 / 1024
             << "MB GPU Memory: " << gpu_memory / 1024 / 1024 << "MB";

    if (total_db_memory > gpu_memory) {
        fold_factor = ceilf((float) total_db_memory / (float) gpu_memory);
    }

    if (gpu_bitcount > 0) {
        unsigned int arg_fold_factor = max_fp_bitcount / gpu_bitcount;
        if (arg_fold_factor < fold_factor) {
            throw std::invalid_argument(
                "GPU bitset not sufficiently small to fit on GPU");
        }
        fold_factor = arg_fold_factor;
    }

    qInfo() << "Putting graphics card data up.";
    if (fold_factor > 1) {
        qDebug() << "Folding databases by at least" << fold_factor
                 << "to fit in gpu memory";
    }

    if (usingGPU()) {
        for (auto db : m_databases) {
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

void GPUSimServer::extractData(const QString& database_fname, int& fp_bitcount,
                               int& fp_count, QString& dbkey,
                               vector<vector<char>>& fingerprint_data,
                               vector<char*>& smiles_vector,
                               vector<char*>& ids_vector)
{
    QFile file(database_fname);
    file.open(QIODevice::ReadOnly);
    QDataStream datastream(&file);
    // Set version so that files will be usable cross-release
    datastream.setVersion(QDataStream::Qt_5_2);
    int version;
    datastream >> version;
    if (version != DATABASE_VERSION) {
        throw std::runtime_error(
            "Database version incompatible with this GPUSim version");
    }

    char* dbkey_char;
    datastream >> dbkey_char;
    dbkey = dbkey_char;
    delete[] dbkey_char;
    datastream >> fp_bitcount;
    datastream >> fp_count;

    int fp_qba_count;
    datastream >> fp_qba_count;
    fingerprint_data.resize(fp_qba_count);
    int current_qba = 1;
    QThreadPool thread_pool;

    for (auto& fp_vector : fingerprint_data) {
        qDebug() << "  loading FP " << current_qba++ << "of" << fp_qba_count;
        auto thread = new DecompressAssignFPRunnable(fp_vector);
        datastream >> thread->compressed_data;
        thread_pool.start(thread);
    }

    int smi_qba_count;
    datastream >> smi_qba_count;
    vector<vector<char*>> smiles_data;
    smiles_data.resize(smi_qba_count);
    current_qba = 1;
    for (auto& local_vector : smiles_data) {
        qDebug() << "  loading SMI " << current_qba++ << "of" << smi_qba_count;
        auto thread = new DecompressAssignStringRunnable(local_vector);
        datastream >> thread->compressed_data;
        thread_pool.start(thread);
    }

    int id_qba_count;
    datastream >> id_qba_count;
    vector<vector<char*>> ids_data;
    ids_data.resize(id_qba_count);
    current_qba = 1;
    for (auto& local_vector : ids_data) {
        qDebug() << "  loading ID " << current_qba++ << "of" << id_qba_count;
        auto thread = new DecompressAssignStringRunnable(local_vector);
        datastream >> thread->compressed_data;
        thread_pool.start(thread);
    }

    qDebug() << "  waiting for data processing threads to finish...";
    thread_pool.waitForDone();

    qDebug() << "  merging smiles vectors";
    for (auto& local_vector : smiles_data) {
        smiles_vector.insert(smiles_vector.end(), local_vector.begin(),
                             local_vector.end());
        local_vector.clear();
    }

    qDebug() << "  merging ID vectors";
    for (auto& local_vector : ids_data) {
        ids_vector.insert(ids_vector.end(), local_vector.begin(),
                          local_vector.end());
        local_vector.clear();
    }

    qDebug() << "  finished merging vectors";
}

bool GPUSimServer::setupSocket()
{
    const QString socket_name("gpusimilarity");
    auto server = new QLocalServer(this);
    if (!server->listen(socket_name)) {
        QString socket_location = QString("/tmp/%1").arg(socket_name);
        QFile::remove(socket_location);
        if (!server->listen(socket_name)) {
            qDebug() << "Server start failed on" << socket_location;
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
                                    const QString& dbname, const QString& dbkey,
                                    unsigned int max_return_count,
                                    float similarity_cutoff,
                                    CalcType calc_type,
                                    vector<char*>& results_smiles,
                                    vector<char*>& results_ids,
                                    vector<float>& results_scores,
                                    unsigned long& approximate_result_count)
{
    if (calc_type == CalcType::GPU) {
        m_databases[dbname]->search(reference, dbkey, max_return_count,
                similarity_cutoff, results_smiles, results_ids,
                results_scores, approximate_result_count);
    } else {
        m_databases[dbname]->search_cpu(reference, dbkey, max_return_count,
                similarity_cutoff, results_smiles, results_ids,
                results_scores, approximate_result_count);
    }
};

void GPUSimServer::newConnection()
{
    auto server = dynamic_cast<QLocalServer*>(sender());
    QLocalSocket* clientConnection = server->nextPendingConnection();

    QObject::connect(clientConnection, &QLocalSocket::disconnected,
                     clientConnection, &QLocalSocket::deleteLater);
    QObject::connect(clientConnection, &QLocalSocket::readyRead, this,
                     &GPUSimServer::incomingSearchRequest);
}

void GPUSimServer::searchDatabases(
        const Fingerprint& query, int results_requested, float similarity_cutoff,
        map<QString, QString>& dbname_to_key, vector<char*>& results_smiles,
        vector<char*>& results_ids, vector<float>& results_scores,
        unsigned long& approximate_result_count)
{
    typedef pair<char*, char*> ResultData;
    typedef pair<float, ResultData> SortableResult;

    vector<SortableResult> sortable_results;
    for (auto name_key_pair : dbname_to_key) {
        const auto& local_dbname = name_key_pair.first;
        const auto& local_key = name_key_pair.second;

        vector<char*> l_results_smiles, l_results_ids;
        vector<float> l_results_scores;
        if (!m_databases.contains(local_dbname)) {
            qDebug() << "Unknown database " << local_dbname << " requested.";
            continue;
        }
        unsigned long local_approximate_result_count;
        similaritySearch(query, local_dbname, local_key,  results_requested,
                         similarity_cutoff,
                         usingGPU() ? CalcType::GPU : CalcType::CPU,
                         l_results_smiles, l_results_ids, l_results_scores,
                         local_approximate_result_count);
        approximate_result_count += local_approximate_result_count;
        for (unsigned int i = 0; i < l_results_smiles.size(); i++) {
            sortable_results.push_back(SortableResult(
                l_results_scores[i],
                ResultData(l_results_smiles[i], l_results_ids[i])));
        }
    }
    std::sort(sortable_results.begin(), sortable_results.end());
    std::reverse(sortable_results.begin(), sortable_results.end());

    map<string, string> smiles_to_ids;
    for (auto result : sortable_results) {
        const auto& smiles = result.second.first;
        const auto& id = result.second.second;
        if (smiles_to_ids.count(smiles) > 0) {
            std::stringstream ss;
            ss << smiles_to_ids[smiles];
            ss << ";:;";
            ss << id;
            smiles_to_ids[smiles] = ss.str();
        } else {
            smiles_to_ids[smiles] = id;
        }
        if (smiles_to_ids.size() >= (unsigned int) results_requested)
            break;
    }

    int written_count = 0;
    set<string> smiles_written;
    for (auto result : sortable_results) {
        const auto& score = result.first;
        const auto& smiles = result.second.first;
        if (smiles_written.count(smiles) > 0)
            continue;
        smiles_written.insert(smiles);
        results_scores.push_back(score);
        results_smiles.push_back(smiles);
        // strdup to hand memory handling off to receiver
        results_ids.push_back(strdup(smiles_to_ids[smiles].c_str()));
        if (++written_count >= results_requested)
            break;
    }
}

void GPUSimServer::incomingSearchRequest()
{
    auto clientConnection = static_cast<QLocalSocket*>(sender());

    // Read incoming Fingerprint binary and put it in Fingerprint object
    QByteArray data = clientConnection->readAll();
    QDataStream qds(&data, QIODevice::ReadOnly);

    int database_search_count;
    qds >> database_search_count;
    map<QString, QString> dbname_to_key;
    for (int i = 0; i < database_search_count; i++) {
        char* dbname;
        char* dbkey;
        qds >> dbname;
        qds >> dbkey;
        dbname_to_key[QString(dbname)] = QString(dbkey);
        delete[] dbname;
        delete[] dbkey;
    }

    int results_requested;
    // Used to guarantee Python reading the pipe reads correct request
    int request_num;
    qds >> request_num;

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
    vector<char*> results_smiles, results_ids;
    vector<float> results_scores;

    QTime timer;
    timer.start();

    unsigned long approximate_result_count = 0;
    searchDatabases(query, results_requested, similarity_cutoff, dbname_to_key,
            results_smiles, results_ids, results_scores,
            approximate_result_count);

    qDebug() << "Search completed, time elapsed:"
             << (float) timer.elapsed() / 1000.0f;

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
    QByteArray ints_qba;
    QDataStream ints_qds(&ints_qba, QIODevice::WriteOnly);
    ints_qds << request_num;
    ints_qds << (int) results_smiles.size();
    ints_qds << (quint64) approximate_result_count;

    clientConnection->write(ints_qba);
    clientConnection->write(output_smiles);
    clientConnection->write(output_ids);
    clientConnection->write(output_scores);
    clientConnection->flush();
}

Fingerprint GPUSimServer::getFingerprint(const int index, const QString& dbname)
{
    return m_databases[dbname]->getFingerprint(index);
}

} // namespace gpusim
