/* -------------------------------------------------------------------------
 * Declares gpusim::FingerprintDB CUDA enabled similarity scoring
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

namespace gpusim
{
class FingerprintDBPriv;

class FingerprintDB
{
  public:
    FingerprintDB(int fp_bitcount, int fp_count, const char* data,
                  std::vector<char*>& smiles_vector,
                  std::vector<char*>& ids_vector);

    /**
     * @brief
     * Copy fingerprint memory up to the GPU, folding to a smaller size
     * if necessary
     *
     * @param fold_factor: Minimum factor to fold fingerprints by, might need
     *                     fold by a bigger factor to get even folding
     */
    void copyToGPU(unsigned int fold_factor);

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
     * @param return_count: Maximum number of results to return
     * @param similarity_cutoff: Minimum similarity score to return
     */
    void search(const Fingerprint& query, std::vector<char*>& results_smiles,
                std::vector<char*>& results_ids,
                std::vector<float>& results_scores,
                unsigned int return_count,
                float similarity_cutoff) const;

    void search_cpu(const Fingerprint& query,
            std::vector<char*>& results_smiles,
            std::vector<char*>& results_ids,
            std::vector<float>& results_scores,
            unsigned int return_count,
            float similarity_cutoff) const;


    char* getSmiles(int index) const { return m_smiles[index]; }
    char* getID(int index) const { return m_ids[index]; }

    size_t getFingerprintDataSize() const { return m_data_size; };
    int getFingerprintBitcount() const { return m_fp_intsize * sizeof(int) * 8; }

  protected:
    // INTERNAL:  A CPU implementation of tanimoto similarity for validation
    float tanimoto_similarity_cpu(const Fingerprint& fp1,
                                  const Fingerprint& fp2) const;

    std::vector<int> fold_data(const std::vector<int>& unfolded) const;

    std::vector<int> m_data, m_folded_data;
    FingerprintDBPriv* m_priv; // Used to conceal cuda types
    int m_count, m_fp_intsize, m_fold_factor;
    size_t m_data_size;
    std::vector<char*> m_smiles;
    std::vector<char*> m_ids;

};

size_t get_available_gpu_memory();

/**
 * @brief
 * Used for CPU sorting to just get top results in O(number_required*N)
 */
void top_results_bubble_sort(std::vector<int>& indices,
        std::vector<float>& scores, int number_required);

std::vector<int> fold_fingerprint(std::vector<int> &, const int);

} // namespace gpusim

#endif
