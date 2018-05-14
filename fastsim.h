#ifndef FASTSIMSERVER_H
#define FASTSIMSERVER_H

#include <memory>
#include <QObject>
#include "types.h"

class QByteArray;
class QSize;
class QString;
class QLocalServer;

namespace fastsim
{
class FingerprintDB;
enum class CalcType {GPU, CPU};

class FastSimServer : public QObject
{
  public:
    /**
     * @brief
     * Start a FastSimServer running on a local socket to communicate with
     * clients.  The server will contain the data in the provided .fsim file.
     * @param database_fname: .fsim file storing relevant data in binary format
     */
    FastSimServer(const QString& database_fname);

    /**
     * @brief
     * Finds the <RETURN_COUNT> most similar compounds stored in the database
     * to the reference fingerprint provided, where RETURN_COUNT is a constant
     * defined in types.h
     *
     * @param query: Fingerprint to find closest matches to
     * @param results_smiles: Vector to store smiles of results
     * @param results_ids: Vector to store IDs of results
     * @param results_scores: Vector to store scores of results
     */
    void similaritySearch(const Fingerprint& reference,
                          std::vector<char*>& results_smiles,
                          std::vector<char*>& results_ids,
                          std::vector<float>& results_scores,
                          CalcType calc_type=CalcType::GPU);

    /**
     * @brief
     * Allows you to fetch a fingerprint from the underlying DB,
     * most useful for testing.
     */
    Fingerprint getFingerprint(const int index);

  public slots:
    void newConnection();

    /**
     * @brief
     * Read a fingerprint from the client, run comparison against the DB,
     * serialize and sends results to the client.
     */
    void incomingSearchRequest();

  private:
    std::shared_ptr<FingerprintDB> m_database;
    std::shared_ptr<QLocalServer> m_server;

    bool setupSocket(const QString& database_fname);
    void extractData(const QString& database_fname, QSize& fingerprint_size,
                     QByteArray& fingerprint_data,
                     std::vector<char*>& smiles_vector,
                     std::vector<char*>& ids_vector);
};
}

#endif
