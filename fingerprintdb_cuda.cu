/* -------------------------------------------------------------------------
 * Implements gpusim::FingerprintDB CUDA enabled similarity
 * scoring
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "fingerprintdb_cuda.h"

#include <iostream>
#include <cmath>

#include <algorithm>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

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

    TanimotoFunctor(const DFingerprint& ref_fp, int fp_intsize, const device_vector<int>& dbdata) : m_ref_fp(ref_fp.data().get()),m_fp_intsize(fp_intsize),m_dbdata(dbdata.data().get())
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

        return (float)common / (float)(total-common);
    };
};


/**
 * @internal
 * Functor used to perform tanimoto similarity on CPU via std::transform
 */
class TanimotoFunctorCPU{
    public:

    const int* m_ref_fp;
    const int m_fp_intsize;
    const int* m_dbdata;

    TanimotoFunctorCPU(const Fingerprint& ref_fp, int fp_intsize, const vector<int>& dbdata) : m_ref_fp(ref_fp.data()),m_fp_intsize(fp_intsize),m_dbdata(dbdata.data())
        {};

    float operator()(const int& fp_index) const
    {
        int total = 0;
        int common = 0;
        int offset = m_fp_intsize*fp_index;
        for(int i=0; i<m_fp_intsize; i++) {
            const int fp1 = m_ref_fp[i];
            const int fp2 = m_dbdata[offset+i];
            total += __builtin_popcount(fp1) + __builtin_popcount(fp2); 
            common += __builtin_popcount(fp1 & fp2);
        }

        return (float)common / (float)(total-common);
    };
};


static inline void swap(vector<int>& indices, vector<float>& scores,
        const int idx1, const int idx2)
{
    int temp = indices[idx1];
    indices[idx1] = indices[idx2];
    indices[idx2] = temp;

    float tempf = scores[idx1];
    scores[idx1] = scores[idx2];
    scores[idx2] = tempf;
}


/**
 * @internal
 * This performs a partial bubble sort, concluding after the top N scores
 * have been sorted.
 * NOTE:  Resulting vectors are *UNSORTED* beyond N positions
 * This version of bubble sort is only O(N*len(scores)), where N is small
 */
void top_results_bubble_sort(vector<int>& indices, vector<float>& scores,
        int number_required)
{
    const int count = indices.size();
    for(int i=0; i<number_required; i++) {
        for(int j=0; j<(count-i-1); j++) {
            if(scores[j] > scores[j+1]) {
                swap(indices, scores, j, j+1);
            }
        }
    }
}


class FingerprintDBPriv
{
    public:
        std::vector<int> data;
        device_vector<int> d_data;

};


FingerprintDB::FingerprintDB(int fp_bitcount, int fp_count, const char* data,
            vector<char*>& smiles_vector,
            std::vector<char*>& ids_vector)
{
    m_priv = new FingerprintDBPriv();
    m_fp_intsize = fp_bitcount / (sizeof(int)*8);  //ASSUMES INT-DIVISIBLE SIZE
    m_count = fp_count;

    int data_size = (fp_bitcount/(sizeof(int)*8))*m_count;
    const int* int_data = (const int*)data;
    m_priv->data.assign(int_data, int_data+data_size);
    m_priv->d_data = m_priv->data;

    // Optimization, take the underlying storage of the incoming vectors, 
    // which won't be used again in calling code
    m_smiles.swap(smiles_vector);
    m_ids.swap(ids_vector);

}


Fingerprint FingerprintDB::getFingerprint(unsigned int index) const
{
    Fingerprint output(m_fp_intsize);

    unsigned int offset = index*m_fp_intsize;
    for(int i=0; i<m_fp_intsize; i++) {
        output[i] = m_priv->data[offset+i];
    }

    return output;
}


void FingerprintDB::search (const Fingerprint& query,
        std::vector<char*>& results_smiles,
        std::vector<char*>& results_ids,
        std::vector<float>& results_scores, unsigned int return_count) const
{
    device_vector<int> d_results_indices(count());
    device_vector<float> d_results_scores(count());

    try
    {
    // Fill indices [0->N), which will be sorted along with scores at end
    thrust::sequence(d_results_indices.begin(), d_results_indices.end());

    // Copy the query fingerprint up to the GPU
    DFingerprint d_ref_fp = query;

    // Use Tanimoto to score similarity of all compounds to query fingerprint
    thrust::transform(d_results_indices.begin(), d_results_indices.end(),
            d_results_scores.begin(),
            TanimotoFunctor(d_ref_fp, m_fp_intsize, m_priv->d_data));

    // Sort scores & indices vectors descending on score
    thrust::sort_by_key(d_results_scores.begin(), d_results_scores.end(),
            d_results_indices.begin(), thrust::greater<float>());
    } catch(thrust::system_error e) {
        std::cerr << "Error!  " << e.what() << std::endl;
    }

    // Push top return_count results to CPU results vectors to be returned
    for(unsigned int i=0;i<return_count;i++) {
        results_smiles.push_back(m_smiles[d_results_indices[i]]);
        results_ids.push_back(m_ids[d_results_indices[i]]);
    }
    results_scores.assign(d_results_scores.begin(),
            d_results_scores.begin()+return_count);

}

void FingerprintDB::search_cpu (const Fingerprint& query,
        std::vector<char*>& results_smiles,
        std::vector<char*>& results_ids,
        std::vector<float>& results_scores, unsigned int return_count) const
{
    const int total = count();
    vector<int> indices(total);
    vector<float> scores(total);

    for(int i=0; i<total; i++) {
        indices[i] = i;
    }

    // Use Tanimoto to score similarity of all compounds to query fingerprint
    std::transform(indices.begin(), indices.end(), scores.begin(),
            TanimotoFunctorCPU(query, m_fp_intsize, m_priv->data));

    top_results_bubble_sort(indices, scores, return_count);

    // Push top return_count results to CPU results vectors to be returned
    for(unsigned int i=total-1;i>=total-return_count;i--) {
        results_smiles.push_back(m_smiles[indices[i]]);
        results_ids.push_back(m_ids[indices[i]]);
        results_scores.push_back(scores[i]);
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

std::vector<int> fold_fingerprint(std::vector<int> &fp, const int factor)
{
    vector<int> new_fp(fp.size()/factor);
    const int INT_SIZE = sizeof(int) * 8;
    const int original_size = INT_SIZE * fp.size();
    // Make sure the new_size is always int-sized
    const int new_size = INT_SIZE * (fp.size() / factor);
    // resize here
    for(int pos=0; pos < original_size; pos++) {
        int int_offset = pos / INT_SIZE;
        int inner_pos = pos % INT_SIZE;
        int bit_on = (fp[int_offset] & (0x01 << inner_pos)) ? 1 : 0;

        int new_pos = pos % new_size;
        int new_int_offset = new_pos / INT_SIZE;
        int new_inner_pos = new_pos % INT_SIZE;
        new_fp[new_int_offset] |= (1 << new_inner_pos) * bit_on;
    }

    return new_fp;
}

} // namespace gpusim
