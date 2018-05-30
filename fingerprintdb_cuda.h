/* -------------------------------------------------------------------------
 * Declares fastsim::FingerprintDB CUDA enabled similarity scoring
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifndef FINGERPRINTDB_CUDA
#define FINGERPRINTDB_CUDA

#include <vector>
#include <utility>
#include "types.h"

class QByteArray;
class QSize;

namespace fastsim
{
class FingerprintDBPriv;

class FingerprintDB
{
  public:
    FingerprintDB(int fp_bitcount, int fp_count, const char* data,
                  std::vector<char*>& smiles_vector,
                  std::vector<char*>& ids_vector);

    // Total number of fingerprints in DB
    unsigned int count() const { return m_count; };

    /**
     * @brief
     * Get an indexed fingerprint that's currently stored in the database
     *
     * @param index: index of fingerprint to return
     * @return: Fingerprint data
     */
    Fingerprint getFingerprint(unsigned int index) const;

    /**
     * @brief
     * Search the database for top ten most similar fingerprints to query
     *
     * @param query: Fingerprint to find closest matches to
     * @param results_smiles: Vector to store smiles of results
     * @param results_ids: Vector to store IDs of results
     * @param results_scores: Vector to store scores of results
     */
    void search(const Fingerprint& query, std::vector<char*>& results_smiles,
                std::vector<char*>& results_ids,
                std::vector<float>& results_scores,
                unsigned int return_count) const;

    void search_cpu (const Fingerprint& query,
            std::vector<char*>& results_smiles,
            std::vector<char*>& results_ids,
            std::vector<float>& results_scores,
            unsigned int return_count) const;


    char* getSmiles(int index) { return m_smiles[index]; }
    char* getID(int index) { return m_ids[index]; }

  protected:
    // INTERNAL:  A CPU implementation of tanimoto similarity for validation
    float tanimoto_similarity_cpu(const Fingerprint& fp1,
                                  const Fingerprint& fp2) const;

    FingerprintDBPriv* m_priv;
    int m_count, m_fp_intsize;
    std::vector<char*> m_smiles;
    std::vector<char*> m_ids;

};

/**
 * @brief
 * Used for CPU sorting to just get top results in O(number_required*N)
 */
void top_results_bubble_sort(std::vector<int>& indices,
        std::vector<float>& scores, int number_required);


} // namespace fastsim

#endif
