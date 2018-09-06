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

#include <QtConcurrent/QtConcurrentRun>
#include <QDebug>
#include <QFuture>
#include <QMutex>

using std::make_shared;
using std::pair;
using std::vector;
using thrust::device_vector;

namespace gpusim
{

unsigned int get_gpu_count()
{
    static int device_count = 0;
    static bool initialized = false;

    if(!initialized) cudaGetDeviceCount(&device_count);
    initialized = true;
    return device_count;
}

unsigned int get_next_gpu()
{
    static int next_device = 0;
    return next_device++ % get_gpu_count(); // Divide by 0 if called w/o GPU
}

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
        float score = static_cast<float>(common) / static_cast<float>(total-common);
        return score >= m_similarity_cutoff ? score : 0;
    };
};


class FingerprintDBPriv
{
    public:
        std::shared_ptr<device_vector<int> > d_data;

};

FingerprintDBStorage::FingerprintDBStorage(FingerprintDB* parent, std::vector<char>& fp_data, 
            int index_offset, int fp_bitcount) : m_parent(parent), m_index_offset(index_offset), 
                                                 m_count(fp_data.size() / (fp_bitcount / 8)), m_gpu_device(get_next_gpu())
{
    const int* int_data = (const int*)fp_data.data();
    const size_t int_size = fp_data.size() / sizeof(int);
    m_data.assign(int_data, int_data+int_size);

}


FingerprintDB::FingerprintDB(int fp_bitcount, int fp_count,
            vector<vector<char> >& data,
            vector<char*>& smiles_vector,
            std::vector<char*>& ids_vector)
{

    m_fp_intsize = fp_bitcount / (sizeof(int)*8);  //ASSUMES INT-DIVISIBLE SIZE
    m_total_count = fp_count;

    int current_fp_count = 0;
    for(unsigned int i=0; i<data.size(); i++) {
        auto& dataset = data[i];
        auto storage = make_shared<FingerprintDBStorage>(this, dataset,
                    current_fp_count, fp_bitcount);
        storage->m_priv = make_shared<FingerprintDBPriv>();
        m_storage.push_back(storage);
        current_fp_count += storage->m_data.size() / m_fp_intsize;
    }

    if(current_fp_count != m_total_count) {
        throw std::runtime_error("Mismatch between FP count and data, potential database corruption.");
    }

    m_total_data_size = (size_t)m_total_count*(size_t)m_fp_intsize*sizeof(int);
    qDebug() << "Database loaded with " << m_total_count << " molecules";

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

    if(m_fold_factor == 1) {
        for(auto& storage : m_storage) {
            cudaSetDevice(storage->m_gpu_device);
            // Have to create vector where correct cuda device is set
            storage->m_priv->d_data = make_shared< device_vector<int> >();
            *(storage->m_priv->d_data) = storage->m_data;
        }
    } else {
        for(auto& storage : m_storage) {
            cudaSetDevice(storage->m_gpu_device);
            // Have to create vector where correct cuda device is set
            storage->m_priv->d_data =  make_shared<device_vector<int> >();
            auto folded_data = fold_data(storage->m_data);
            *(storage->m_priv->d_data) = folded_data;
        }
    }
}


Fingerprint FingerprintDB::getFingerprint(unsigned int index) const
{
    Fingerprint output(m_fp_intsize);
    int slice_index_offset=0;
    auto storage = m_storage[0];
    for(unsigned int i=1; i<m_storage.size(); i++) {
        if(m_storage[i]->m_index_offset >= index) break;
        storage = m_storage[i];
        slice_index_offset = storage->m_index_offset;
    }
    index -= slice_index_offset;

    unsigned int offset = index*m_fp_intsize;
    for(int i=0; i<m_fp_intsize; i++) {
        output[i] = storage->m_data[offset+i];
    }

    return output;
}


void FingerprintDB::search_storage(const Fingerprint& query,
        const std::shared_ptr<FingerprintDBStorage>& storage,
        vector<SortableResult>* sortable_results,
        unsigned int return_count,
        float similarity_cutoff) const
{
    cudaSetDevice(storage->m_gpu_device);
    static QMutex mutex;
    vector<int> indices;
    std::vector<char*> results_smiles;
    std::vector<char*> results_ids;
    std::vector<float> results_scores;
    device_vector<float> d_results_scores(storage->m_count);
    device_vector<int> d_results_indices(storage->m_count);
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

        const int folded_fp_intsize = m_fp_intsize / m_fold_factor;
        // Use Tanimoto to score similarity of all compounds to query fingerprint
        thrust::transform(d_results_indices.begin(), d_results_indices.end(),
                d_results_scores.begin(),
                TanimotoFunctor(d_ref_fp, folded_fp_intsize, *(storage->m_priv->d_data),
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

        int results_to_consider = 0;
        results_to_consider = std::min(indices_size, return_count*m_fold_factor);

        indices.assign(d_results_indices.begin(), 
                d_results_indices.begin()+results_to_consider);

    } catch(thrust::system_error e) {
        qDebug() << "Error!  " << e.what();
    }

    if(m_fold_factor == 1) { // If we don't fold, we can take exact GPU results
        // Push top return_count results to CPU results vectors to be returned
        for(auto index : indices) {
            int offset_index = index + storage->m_index_offset;
            results_smiles.push_back(m_smiles[offset_index]);
            results_ids.push_back(m_ids[offset_index]);
        }
        results_scores.assign(d_results_scores.begin(),
                d_results_scores.begin()+indices.size());
    } else { // If we folded, we need to recalculate scores with full fingerprints
        results_scores.resize(indices.size());
        for(unsigned int i=0;i<indices.size();i++) {
            int offset_index = indices[i] + storage->m_index_offset;
            results_scores[i] = tanimoto_similarity_cpu(query,
                    getFingerprint(offset_index));
            // Uncomment below to debug pre vs post folding scores
            // qDebug() << results_scores[i] << " vs " << d_results_scores[i];
        }
        top_results_bubble_sort(indices, results_scores, return_count);

        return_count = std::min((size_t)return_count, indices.size());
        results_scores.resize(return_count);
        for(unsigned int i=0;i<return_count;i++) {
            // Check whether the re-scored similarity is too low
            if(results_scores[i] < similarity_cutoff) {
                results_scores.resize(i);
                break;
            }
            results_ids.push_back(m_ids[indices[i]+storage->m_index_offset]);
            results_smiles.push_back(m_smiles[indices[i]+storage->m_index_offset]);
        }
    }

    mutex.lock();
    for(unsigned int i=0; i<results_smiles.size(); i++) {
        sortable_results->push_back(SortableResult(results_scores[i], 
                    ResultData(results_smiles[i], results_ids[i])));
    }
    mutex.unlock();
}


void FingerprintDB::search(const Fingerprint& query,
        std::vector<char*>& results_smiles,
        std::vector<char*>& results_ids,
        std::vector<float>& results_scores,
        unsigned int return_count,
        float similarity_cutoff) const
{
    vector<SortableResult> sortable_results;

    vector<QFuture<void> > futures;
    for(auto& storage : m_storage) {
        QFuture<void> future = QtConcurrent::run(this,
                &FingerprintDB::search_storage, query, storage,
                &sortable_results, return_count, similarity_cutoff);
        futures.push_back(future);
    }
    for(auto& future : futures) {
        future.waitForFinished();
    }
    std::sort(sortable_results.begin(), sortable_results.end());
    std::reverse(sortable_results.begin(), sortable_results.end());

    for(auto result : sortable_results) {
        results_scores.push_back(result.first);
        results_smiles.push_back(result.second.first);
        results_ids.push_back(result.second.second);
    }
    int result_size = std::min((int)return_count, (int)results_scores.size());
    results_scores.resize(result_size);
    results_smiles.resize(result_size);
    results_ids.resize(result_size);
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
    for(unsigned int i=0; i<get_gpu_count(); i++) {
        size_t lfree;
        cudaSetDevice(i);
        cudaMemGetInfo(&lfree, &total);
        free += lfree;
    }

    // Comment out below line to force-test folding:
    // free = 100*1024*1024;

    return free;
}

} // namespace gpusim
