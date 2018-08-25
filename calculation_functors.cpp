#include "calculation_functors.h"


namespace gpusim
{

void TanimotoFunctorCPU::operator()(const int& fp_index) const
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

    m_output[fp_index] = static_cast<float>(common) / static_cast<float>(total-common);
}

void FoldFingerprintFunctorCPU::operator()(const int& fp_index) const
{
    const int* fp =  &(m_unfolded[fp_index * m_unfolded_fp_intsize]);
    int* new_fp = &(m_folded[fp_index * m_folded_fp_intsize]);
    const int INT_SIZE = sizeof(int) * 8;
    const int original_size = INT_SIZE * m_unfolded_fp_intsize;
    // Make sure the new_size is always int-sized
    const int new_size = INT_SIZE * m_folded_fp_intsize;
    // resize here
    for(int pos=0; pos < original_size; pos++) {
        int int_offset = pos / INT_SIZE;
        int inner_pos = pos % INT_SIZE;
        int bit_on = (fp[int_offset] & (0x01 << inner_pos)) ? 1 : 0;

        int new_pos = pos % new_size;
        int new_int_offset = new_pos / INT_SIZE;
        const int new_inner_pos = inner_pos; // Always the same
        new_fp[new_int_offset] |= (1 << new_inner_pos) * bit_on;
    }

}

}
