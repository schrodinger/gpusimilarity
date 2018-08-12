/* -------------------------------------------------------------------------
 * Implements gpusim::FingerprintDB CUDA enabled similarity
 * scoring
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include <cuda_runtime.h>
#include "fingerprintdb_cuda.h"

#include <iostream>
#include <cmath>

#include <algorithm>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

using std::cerr;
using std::endl;
using std::vector;
using thrust::device_vector;


namespace gpusim
{

typedef device_vector<int> DFingerprint;

/**
 * @internal
 * Functor used to perform tanimoto similarity on GPGPU via thrust::transform
 */
struct TanimotoFunctor {

    const int* m_ref_fp;
    const int m_fp_intsize;
    const int* m_dbdata;
    const float m_similarity_cutoff;

    TanimotoFunctor(const DFingerprint& ref_fp, int fp_intsize,
            const device_vector<int>& dbdata, float similarity_cutoff) :
        m_ref_fp(ref_fp.data().get()),m_fp_intsize(fp_intsize),m_dbdata(dbdata.data().get()),
        m_similarity_cutoff(similarity_cutoff)
        {};

    __device__ float
    operator()(const int& fp_index) const
    {
        int total = 0;
        int common = 0;
        int offset = m_fp_intsize*fp_index;
        for(int i=0; i<m_fp_intsize; i++) {
            const int fp1 = m_ref_fp[i];
            const int fp2 = m_dbdata[offset+i];
            total += __popc(fp1) + __popc(fp2); 
            common += __popc(fp1 & fp2);
        }
        float score = (float)common / (float)(total-common);
        return score > m_similarity_cutoff ? score : 0;
    };
};


class FingerprintDBPriv
{
    public:
        device_vector<int> d_data;

};


FingerprintDB::FingerprintDB(int fp_bitcount, int fp_count, const char* data,
            vector<char*>& smiles_vector,
            std::vector<char*>& ids_vector)
{

    m_priv = new FingerprintDBPriv();
    m_fp_intsize = fp_bitcount / (sizeof(int)*8);  //ASSUMES INT-DIVISIBLE SIZE
    m_count = fp_count;

    m_data_size = m_fp_intsize*m_count;
    const int* int_data = (const int*)data;
    m_data.assign(int_data, int_data+m_data_size);

    // Optimization, take the underlying storage of the incoming vectors, 
    // which won't be used again in calling code
    m_smiles.swap(smiles_vector);
    m_ids.swap(ids_vector);

}


void FingerprintDB::copyToGPU(unsigned int fold_factor)
{
    m_fold_factor = fold_factor;
    while(m_fp_intsize % m_fold_factor != 0) {
        m_fold_factor++;
    }

    if(fold_factor == 1) {
        m_priv->d_data = m_data;
    } else {
        m_folded_data = fold_data(m_data);
        m_priv->d_data = m_folded_data;
    }
}


Fingerprint FingerprintDB::getFingerprint(unsigned int index) const
{
    Fingerprint output(m_fp_intsize);

    unsigned int offset = index*m_fp_intsize;
    for(int i=0; i<m_fp_intsize; i++) {
        output[i] = m_data[offset+i];
    }

    return output;
}


void FingerprintDB::search(const Fingerprint& query,
        std::vector<char*>& results_smiles,
        std::vector<char*>& results_ids,
        std::vector<float>& results_scores,
        unsigned int return_count,
        float similarity_cutoff) const
{
    device_vector<int> d_results_indices(count());
    device_vector<float> d_results_scores(count());
    vector<int> indices(return_count*m_fold_factor);
    int results_to_consider = 0;

    try
    {
    // Fill indices [0->N), which will be sorted along with scores at end
    thrust::sequence(d_results_indices.begin(), d_results_indices.end());

    DFingerprint d_ref_fp;
    if(m_fold_factor == 1) {
        // Copy the query fingerprint up to the GPU
        d_ref_fp = query;
    } else {
        auto folded = fold_data(query);
        d_ref_fp = folded;
    }

    int folded_fp_intsize = m_fp_intsize / m_fold_factor;
    // Use Tanimoto to score similarity of all compounds to query fingerprint
    thrust::transform(d_results_indices.begin(), d_results_indices.end(),
            d_results_scores.begin(),
            TanimotoFunctor(d_ref_fp, folded_fp_intsize, m_priv->d_data,
                similarity_cutoff));

    auto indices_end = d_results_indices.end();
    auto scores_end = d_results_scores.end();
    if(similarity_cutoff > 0) {
        indices_end = thrust::remove_if(d_results_indices.begin(),
                d_results_indices.end(), d_results_scores.begin(),
                thrust::logical_not<bool>());
        scores_end = thrust::remove(d_results_scores.begin(),
                d_results_scores.end(), 0);
    }
    unsigned int indices_size = std::distance(d_results_indices.begin(),
            indices_end);

    // Sort scores & indices vectors descending on score
    thrust::sort_by_key(d_results_scores.begin(), scores_end,
            d_results_indices.begin(), thrust::greater<float>());

    results_to_consider = std::min(indices_size, return_count*m_fold_factor);

    indices.assign(d_results_indices.begin(), 
            d_results_indices.begin()+results_to_consider);

    } catch(thrust::system_error e) {
        cerr << "Error!  " << e.what() << endl;
    }


    if(m_fold_factor == 1) { // If we don't fold, we can take exact GPU results
        // Push top return_count results to CPU results vectors to be returned
        for(unsigned int i=0;i<results_to_consider;i++) {
            results_smiles.push_back(m_smiles[d_results_indices[i]]);
            results_ids.push_back(m_ids[d_results_indices[i]]);
        }
        results_scores.assign(d_results_scores.begin(),
                d_results_scores.begin()+results_to_consider);
    } else { // If we folded, we need to recalculate scores with full fingerprints
        results_scores.resize(results_to_consider);
        for(unsigned int i=0;i<results_to_consider;i++) {
            results_scores[i] = tanimoto_similarity_cpu(query,
                    getFingerprint(indices[i]));
            // Uncomment below to debug pre vs post folding scores
            // cerr << results_scores[i] << " vs " << d_results_scores[i] << endl;
        }
        top_results_bubble_sort(indices, results_scores, return_count);

        return_count = std::min((size_t)return_count, indices.size());
        results_scores.resize(return_count);
        for(unsigned int i=0;i<return_count;i++) {
            results_ids.push_back(m_ids[indices[i]]);
            results_smiles.push_back(m_smiles[indices[i]]);
        }
    }
}

/**
 * @brief
 * A CPU implementation of tanimoto similarity, meant purely for testing.
 */
float FingerprintDB::tanimoto_similarity_cpu(const Fingerprint& fp1,
        const Fingerprint& fp2) const
{

    int total = 0;
    int common = 0;
    for(int i=0; i<m_fp_intsize; i++) {
        total += __builtin_popcount(fp1[i]) + __builtin_popcount(fp2[i]); 
        common += __builtin_popcount(fp1[i] & fp2[i]);
    }

    return (float)common / (float)(total-common);
}

size_t get_available_gpu_memory()
{
    size_t free=0, total=0;
    cudaMemGetInfo(&free, &total);
    // Comment out below line to force-test folding:
    // free = 100*1024*1024;

    return free;
}

} // namespace gpusim
