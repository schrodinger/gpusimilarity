/* -------------------------------------------------------------------------
 * Implements gpusim::FingerprintDB *NON*-CUDA functionality
 * This is in its own file so that cuda and Qt don't have to intermix
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "fingerprintdb_cuda.h"

#include <QDir>
#include <QDebug>
#include <QString>
#include <QDataStream>
#include <QThreadPool>
#include <QCryptographicHash>
#include <QtConcurrent/QtConcurrentMap>

#include "calculation_functors.h"

using namespace gpusim;
using std::vector;

void FingerprintDB::search_cpu(const Fingerprint& query, const QString& dbkey,
                               unsigned int max_return_count,
                               float similarity_cutoff,
                               std::vector<char*>& results_smiles,
                               std::vector<char*>& results_ids,
                               std::vector<float>& results_scores,
                               unsigned long& approximate_result_count) const
{
    if (dbkey != m_dbkey) {
        qDebug() << "Key check failed, returning empty results";
        return;
    }
    const int total = count();
    vector<int> indices(total);
    std::iota(indices.begin(), indices.end(), 0);

    vector<float> scores(total);

    // TODO:  CPU searching broken for >1GB, needs to be fixed below
    // TODO:  CPU searching broken, not giving approximate total count back

    // Scoring parallelizes well, but bottleneck is now sorting
    QtConcurrent::blockingMap(
        indices,
        TanimotoFunctorCPU(query, m_fp_intsize, m_storage[0]->m_data, scores));

    top_results_bubble_sort(indices, scores, max_return_count);

    // Push top return_count results to CPU results vectors to be returned
    for (unsigned int i = 0; i < max_return_count; i++) {
        results_smiles.push_back(m_smiles[indices[i]]);
        results_ids.push_back(m_ids[indices[i]]);
        results_scores.push_back(scores[i]);
    }
}

std::vector<int>
FingerprintDB::fold_data(const std::vector<int>& unfolded) const
{
    vector<int> folded(unfolded.size() / m_fold_factor);
    std::fill(folded.begin(), folded.end(), 0);

    vector<int> indices(unfolded.size() / m_fp_intsize);
    std::iota(indices.begin(), indices.end(), 0);

    QtConcurrent::blockingMap(
        indices, FoldFingerprintFunctorCPU(m_fold_factor, m_fp_intsize,
                                           unfolded, folded));
    return folded;
}

QByteArray FingerprintDB::getHash() const
{
    QCryptographicHash algo(QCryptographicHash::Sha256);
    for (const auto& storage: m_storage) {
        const vector<int>& fps = storage->m_data;
        algo.addData(reinterpret_cast<const char*>(fps.data()),
                     sizeof(int)*fps.size());
    }
    return algo.result();
}

QString FingerprintDB::getFFPCacheFilename(unsigned int fold_factor) const
{
    return QString::number(fold_factor) + '-' + getHash().toHex();
}

std::unique_ptr<QFile>
FingerprintDB::openFFPCacheFile(const QString& cache_directory,
                                unsigned int fold_factor) const
{
    std::unique_ptr<QFile> outcome;

    if (cache_directory.isEmpty()) {
        return outcome;
    }

    QDir dir(cache_directory);
    if (!dir.exists()) {
        if (!dir.mkpath(".")) {
            qWarning() << Q_FUNC_INFO << ": could not create"
                       << qUtf8Printable(cache_directory);
            return outcome;
        }
    }

    outcome.reset(new QFile(dir.filePath(getFFPCacheFilename(fold_factor))));

    if (outcome->exists()) {
        if (!outcome->open(QIODevice::ReadOnly)) {
            qWarning() << Q_FUNC_INFO << ": could not open"
                       << qUtf8Printable(outcome->fileName());
            outcome.reset();
        }
    } else {
        if (!outcome->open(QIODevice::WriteOnly)) {
            qWarning() << Q_FUNC_INFO << ": could not create"
                       << qUtf8Printable(outcome->fileName());
            outcome.reset();
        }
    }

    return outcome;
}

void FingerprintDB::copyToGPU
    (unsigned int fold_factor, const QString& cache_directory)
{
    m_fold_factor = fold_factor;
    while (m_fp_intsize % m_fold_factor != 0) {
        m_fold_factor++;
    }

    if (m_fold_factor == 1) {
        for (const auto& storage : m_storage) {
            storage->copyToGPU(storage->m_data);
        }
    } else {
        std::unique_ptr<QFile> cache_file =
            openFFPCacheFile(cache_directory, fold_factor);

        QDataStream cache_stream;
        if (cache_file && cache_file->openMode() != QIODevice::NotOpen) {
            cache_stream.setDevice(cache_file.get());
            cache_stream.setVersion(QDataStream::Qt_5_2);
            qDebug() << Q_FUNC_INFO << ": cache:"
                     << qUtf8Printable(cache_file->fileName())
                     << " mode:" << cache_file->openMode();
        }

        std::vector<int> folded_data;
        for (const auto& storage : m_storage) {
            const size_t folded_size = storage->m_data.size()/m_fold_factor;

            if (cache_file && cache_file->openMode() == QIODevice::ReadOnly) {
                folded_data.resize(folded_size);
                cache_stream.readRawData(
                    (char*) folded_data.data(), folded_size*sizeof(int));
            } else {
                folded_data = fold_data(storage->m_data);
                if (cache_file &&
                    cache_file->openMode() == QIODevice::WriteOnly) {
                    cache_stream.writeRawData(
                        (char*) folded_data.data(), folded_size*sizeof(int));
                }
            }
            storage->copyToGPU(folded_data);
        }
    }
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
    for (int i = 0; i < number_required; i++) {
        for (int j = (count - 1); j > i; j--) {
            if (scores[j] > scores[j - 1]) {
                swap(indices, scores, j, j - 1);
            }
        }
    }
}
} // End namespace gpusim
