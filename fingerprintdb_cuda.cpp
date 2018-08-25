/* -------------------------------------------------------------------------
 * Implements gpusim::FingerprintDB *NON*-CUDA functionality
 * This is in its own file so that cuda and Qt don't have to intermix
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "fingerprintdb_cuda.h"

#include <QThreadPool>
#include <QtConcurrent/QtConcurrentMap>

#include "calculation_functors.h"

using namespace gpusim;
using std::vector;

void FingerprintDB::search_cpu (const Fingerprint& query,
        std::vector<char*>& results_smiles,
        std::vector<char*>& results_ids,
        std::vector<float>& results_scores,
        unsigned int return_count,
        float similarity_cutoff) const
{
    const int total = count();
    vector<int> indices(total);
    std::iota(indices.begin(), indices.end(), 0);

    vector<float> scores(total);

    // Scoring parallelizes well, but bottleneck is now sorting
    QtConcurrent::blockingMap(indices,
            TanimotoFunctorCPU(query, m_fp_intsize, m_data, scores));

    top_results_bubble_sort(indices, scores, return_count);

    // Push top return_count results to CPU results vectors to be returned
    for(unsigned int i=0;i<return_count;i++) {
        results_smiles.push_back(m_smiles[indices[i]]);
        results_ids.push_back(m_ids[indices[i]]);
        results_scores.push_back(scores[i]);
    }
}


std::vector<int> FingerprintDB::fold_data(const std::vector<int>& unfolded) const
{
    vector<int> folded(unfolded.size() / m_fold_factor);
    std::fill(folded.begin(), folded.end(), 0);

    vector<int> indices(unfolded.size() / m_fp_intsize);
    std::iota(indices.begin(), indices.end(), 0);

    QtConcurrent::blockingMap(indices,
            FoldFingerprintFunctorCPU(m_fold_factor, m_fp_intsize, unfolded,
                folded));
    return folded;
}


namespace gpusim
{
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
        for(int j=(count-1); j>i; j--) {
            if(scores[j] > scores[j-1]) {
                swap(indices, scores, j, j-1);
            }
        }
    }
}
} // End namespace gpusim
