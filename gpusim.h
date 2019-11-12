#pragma once

#include <map>
#include <memory>
#include <string>

#include "qstring_hash.h"
#include "types.h"
#include <QHash>
#include <QObject>
#include <QString>

class QByteArray;
class QSize;
class QString;
class QLocalServer;

namespace gpusim
{
class FingerprintDB;
enum class CalcType { GPU, CPU };

class GPUSimServer : public QObject
{
  public:
    /**
     * @brief
     * Start a GPUSimServer running on a local socket to communicate with
     * clients.  The server will contain the data in the provided .fsim file.
     * @param database_fname: .fsim file storing relevant data in binary format
     * @param gpu_bitcount: If 0, then auto-calculate optimal value
     * @param cache_directory: Path to cache directory. Empty string disables
     *     folded fingerprints cacheing.
     */
    explicit GPUSimServer(const QStringList& database_fnames,
                          int gpu_bitcount = 0,
                          const QString& cache_directory = QString());

    /**
     * @brief
     * Finds the most similar compounds stored in the database
     * to the reference fingerprint provided
     *
     * @param reference: Fingerprint to find closest matches to
     * @param dbname: Which database to search against
     * @param dbkey: Key to access that database
     * @param max_return_count: Maximum number of results to return
     * @param similarity_cutoff: Minimum similarity score to return molecules
     * @param calc_type: Whether to search on CPU or GPU
     * @param results_smiles: Vector to store smiles of results
     * @param results_ids: Vector to store IDs of results
     * @param results_scores: Vector to store scores of results for
     */
    void similaritySearch(const Fingerprint& reference, const QString& dbname,
                          const QString& dbkey, unsigned int max_return_count,
                          float similarity_cutoff, CalcType calc_type,
                          std::vector<char*>& results_smiles,
                          std::vector<char*>& results_ids,
                          std::vector<float>& results_scores,
                          unsigned long& approximate_result_count);

    void searchDatabases(const Fingerprint& reference, int results_requested,
                         float similarity_cutoff,
                         std::map<QString, QString>& dbname_to_key,
                         std::vector<char*>& results_smiles,
                         std::vector<char*>& results_ids,
                         std::vector<float>& results_scores,
                         unsigned long& approximate_result_count);

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
    bool usingGPU();

  private:
    QHash<QString, std::shared_ptr<FingerprintDB>> m_databases;
    bool m_use_gpu = true;

    bool setupSocket();
    static void extractData(
        const QString& database_fname,
        int& fp_bitcount,
        int& fp_count,
        QString& dbkey,
        std::vector<std::vector<char>>& fingerprint_data,
        std::vector<char*>& smiles_vector,
        std::vector<char*>& ids_vector);
};
} // namespace gpusim
