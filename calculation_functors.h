#pragma once

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include "cuda_types.h"
#endif
#include <QVector>
#include "types.h"

namespace gpusim
{
/**
 * Functor used to perform tanimoto similarity on GPGPU via thrust::transform
 */
#ifdef __CUDACC__
class TanimotoFunctor {

    const int* m_ref_fp;
    const int m_fp_intsize;
    const int* m_dbdata;

public:

    TanimotoFunctor(const gpusim::DFingerprint& ref_fp, int fp_intsize, const thrust::device_vector<int>& dbdata) : m_ref_fp(ref_fp.data().get()),m_fp_intsize(fp_intsize),m_dbdata(dbdata.data().get())
        {};

 /*
  * This function is defined in calculation_functors.cu as cuda prefers device
  * code defined in location of use.
  */
    __device__ float operator()(const int& fp_index) const;
};
#endif


/**
 * Functor used to perform tanimoto similarity on CPU via std::transform
 */
struct TanimotoFunctorCPU{
    const int* m_ref_fp;
    const int m_fp_intsize;
    const int* m_dbdata;
    float* m_output;

    TanimotoFunctorCPU(const gpusim::Fingerprint& ref_fp, int fp_intsize, const std::vector<int>& dbdata, std::vector<float>& output) : m_ref_fp(ref_fp.data()),m_fp_intsize(fp_intsize),m_dbdata(dbdata.data()),m_output(output.data())
        {};

    void operator()(const int& fp_index) const;
};

/**
 * Functor used to fold fingerprints down on CPU
 */
class FoldFingerprintFunctorCPU{
    const int m_unfolded_fp_intsize;
    const int m_folded_fp_intsize;
    const int* m_unfolded;
    int* m_folded;

public:
    FoldFingerprintFunctorCPU(const int factor, const int fp_intsize,
            const std::vector<int>& unfolded, std::vector<int>& folded) :
        m_unfolded_fp_intsize(fp_intsize),
        m_folded_fp_intsize(fp_intsize/factor), m_unfolded(unfolded.data()),
        m_folded(folded.data())
    {};

    void operator()(const int& fp_index) const;
};

}; // End namespace gpusim
