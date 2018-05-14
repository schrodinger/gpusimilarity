/* -------------------------------------------------------------------------
 * Implements fastsim::FastSimServer
 * Reads a binary .fsim file and creates a FingerprintDB on GPU.  Uses a local
 * socket to communicate with querying processes.  Fingerprint agnostic.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "fastsim.h"

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

using std::shared_ptr;
using std::vector;
using fastsim::FingerprintDB;
using fastsim::Fingerprint;

namespace fastsim
{
FastSimServer::FastSimServer(const QString& database_fname)
{

    // Read from .fsim file into byte arrays
    QSize fingerprint_size;
    QByteArray fingerprint_data;
    vector<char*> smiles_vector;
    vector<char*> ids_vector;
    extractData(database_fname, fingerprint_size, fingerprint_data,
                smiles_vector, ids_vector);

    // Create new FingerprintDB for querying on GPU
    FingerprintDB* fpdb =
        new FingerprintDB(fingerprint_size.width(), fingerprint_size.height(),
                          fingerprint_data.data(), smiles_vector, ids_vector);
    m_database.reset(fpdb);

    if (!setupSocket(database_fname))
        return;
    qInfo() << "Ready for searches.";
};

void FastSimServer::extractData(const QString& database_fname,
                                QSize& fingerprint_size,
                                QByteArray& fingerprint_data,
                                vector<char*>& smiles_vector,
                                vector<char*>& ids_vector)
{
    QFile file(database_fname);
    file.open(QIODevice::ReadOnly);
    QDataStream datastream(&file);
    // Set version so that files will be usable cross-release
    datastream.setVersion(QDataStream::Qt_5_9);

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

bool FastSimServer::setupSocket(const QString& database_fname)
{
    m_server.reset(new QLocalServer());
    QFileInfo file_info(database_fname);
    QString socket_name = file_info.baseName();
    if (!m_server->listen(socket_name)) {
        QString socket_location = QString("/tmp/%1").arg(socket_name);
        QFile::remove(socket_location);
        if (!m_server->listen(socket_name)) {
            qDebug() << "Server start failed on " << socket_location;
            auto app = QCoreApplication::instance();
            app->exit(1);
            return false;
        }
    }

    QLocalServer::connect(m_server.get(), &QLocalServer::newConnection, this,
                          &FastSimServer::newConnection);

    return true;
}

void FastSimServer::similaritySearch(const Fingerprint& reference,
        vector<char*>& results_smiles, vector<char*>& results_ids,
        vector<float>& results_scores, CalcType calc_type)
{
    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, nullptr);

    if(calc_type == CalcType::GPU) {
        m_database->search(reference, results_smiles, results_ids,
                results_scores);
    } else {
        m_database->search_cpu(reference, results_smiles, results_ids,
                results_scores);
    }

    gettimeofday(&tval_after, nullptr);

    timersub(&tval_after, &tval_before, &tval_result);

    // Using printf here because decimal formatting like this is a PITA w/ cout
    printf("Search completed, time elapsed: %ld.%06ld\n",
           (long int) tval_result.tv_sec, (long int) tval_result.tv_usec);
};

void FastSimServer::newConnection()
{
    QLocalSocket* clientConnection = m_server->nextPendingConnection();
    QObject::connect(clientConnection, &QLocalSocket::disconnected,
                     clientConnection, &QLocalSocket::deleteLater);
    QObject::connect(clientConnection, &QLocalSocket::readyRead, this,
                     &FastSimServer::incomingSearchRequest);
}

void FastSimServer::incomingSearchRequest()
{
    auto clientConnection = static_cast<QLocalSocket*>(sender());

    // Read incoming Fingerprint binary and put it in Fingerprint object
    QByteArray data = clientConnection->readAll();
    const int* raw_data = reinterpret_cast<const int*>(data.constData());
    int fp_int_size = data.size() / sizeof(int);

    Fingerprint query(fp_int_size);
    query.assign(raw_data, raw_data + fp_int_size);

    // Perform similarity search and return results in relevant vectors
    vector<char *> results_smiles, results_ids;
    vector<float> results_scores;
    similaritySearch(query, results_smiles, results_ids, results_scores);

    // Create QByteArrays and QDataStreams to write to corresponding arrays
    QByteArray output_smiles, output_ids, output_scores;
    QDataStream smiles_stream(&output_smiles, QIODevice::WriteOnly);
    QDataStream ids_stream(&output_ids, QIODevice::WriteOnly);
    QDataStream scores_stream(&output_scores, QIODevice::WriteOnly);
    for (int i = 0; i < RETURN_COUNT; i++) {
        smiles_stream << results_smiles[i];
        ids_stream << results_ids[i];
        scores_stream << results_scores[i];
    }

    // Transmit binary data to client and flush the buffered data
    clientConnection->write(output_smiles);
    clientConnection->write(output_ids);
    clientConnection->write(output_scores);
    clientConnection->flush();
}


Fingerprint FastSimServer::getFingerprint(const int index)
{
    return m_database->getFingerprint(index);
}

} // end fastsim
