/* -------------------------------------------------------------------------
 * Declares gpusim::FingerprintDB CUDA enabled similarity scoring
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifndef FINGERPRINTDB_CUDA
#define FINGERPRINTDB_CUDA

#include <memory>
#include <vector>
#include <utility>

#include <QObject> 

#include "types.h"

class QByteArray;
class QSize;

namespace gpusim
{

class FingerprintDB;
class FingerprintDBPriv;

typedef std::pair<char*, char*> ResultData;
typedef std::pair<float, ResultData > SortableResult;

unsigned int get_gpu_count();
unsigned int get_next_gpu();

class FingerprintDBStorage
{
public:
    friend class FingerprintDB;

    FingerprintDBStorage(FingerprintDB* parent, std::vector<char>& fp_data, 
            int index_offset, int fp_bitcount);
    unsigned int getOffsetIndex(unsigned int without_offset);

private:
    FingerprintDB* m_parent;
    std::vector<int> m_data;
    std::shared_ptr<FingerprintDBPriv> m_priv; // Used to conceal cuda types
    const unsigned int m_index_offset;
    const int m_count;
    const int m_gpu_device;
};
    
class FingerprintDB : public QObject
{
  friend class FingerprintDBStorage;

  public:
    FingerprintDB(int fp_bitcount, int fp_count,
            std::vector<std::vector<char> >& data,
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
    unsigned int count() const { return m_total_count; };

    /*
     * @brief
     * This function takes an index in the range of all molecules, and finds
     * which storage block that index is inside, and returns that storage block
     * and the local index inside that block
     */
    void getStorageAndLocalIndex(unsigned int offset_index,
            FingerprintDBStorage** storage, unsigned int* local_index) const;

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

    size_t getFingerprintDataSize() const { return m_total_data_size; };
    int getFingerprintBitcount() const { return m_fp_intsize * sizeof(int) * 8; }

    void search_storage(const Fingerprint& query,
            const std::shared_ptr<FingerprintDBStorage>& storage,
            std::vector<SortableResult>* sortable_results,
            unsigned int return_count, float similarity_cutoff) const;

  protected:
    // INTERNAL:  A CPU implementation of tanimoto similarity for validation
    float tanimoto_similarity_cpu(const Fingerprint& fp1,
                                  const Fingerprint& fp2) const;
    std::vector<int> fold_data(const std::vector<int>& unfolded) const;

    std::vector<std::shared_ptr<FingerprintDBStorage> > m_storage;

    int m_total_count, m_fp_intsize, m_fold_factor;
    size_t m_total_data_size;
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
