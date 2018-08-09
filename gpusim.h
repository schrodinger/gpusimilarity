#pragma once

#include <memory>
#include <QHash>
#include <QObject>
#include <QString>
#include "types.h"

class QByteArray;
class QSize;
class QString;
class QLocalServer;

namespace gpusim
{
class FingerprintDB;
enum class CalcType {GPU, CPU};

class GPUSimServer : public QObject
{
  public:
    /**
     * @brief
     * Start a GPUSimServer running on a local socket to communicate with
     * clients.  The server will contain the data in the provided .fsim file.
     * @param database_fname: .fsim file storing relevant data in binary format
     */
    GPUSimServer(const QStringList& database_fnames);

    /**
     * @brief
     * Finds the <return_count> most similar compounds stored in the database
     * to the reference fingerprint provided
     *
     * @param query: Fingerprint to find closest matches to
     * @param dbname: Which database to search against
     * @param results_smiles: Vector to store smiles of results
     * @param results_ids: Vector to store IDs of results
     * @param results_scores: Vector to store scores of results
     */
    void similaritySearch(const Fingerprint& reference,
                          const QString& dbname,
                          std::vector<char*>& results_smiles,
                          std::vector<char*>& results_ids,
                          std::vector<float>& results_scores,
                          unsigned int return_count,
                          CalcType calc_type=CalcType::GPU);

    /**
     * @brief
     * Allows you to fetch a fingerprint from the underlying DB,
     * most useful for testing.
     */
    Fingerprint getFingerprint(const int index, const QString& dbname);

  public slots:
    void newConnection();

    /**
     * @brief
     * Read a fingerprint from the client, run comparison against the DB,
     * serialize and sends results to the client.
     */
    void incomingSearchRequest();
    void setUseGPU(bool use_gpu) { m_use_gpu = use_gpu; }
    bool usingGPU(){ return m_use_gpu; }

  private:
    QHash<QString, std::shared_ptr<FingerprintDB>> m_databases;
    bool m_use_gpu = true;

    bool setupSocket(const QString& socket_name);
    void extractData(const QString& database_fname, QSize& fingerprint_size,
                     QByteArray& fingerprint_data,
                     std::vector<char*>& smiles_vector,
                     std::vector<char*>& ids_vector);
};
}
