/*
 * NOTE:  I mainly left this in here because this is working with waf, and
 * I didn't want to have to figure out how to do this again if problems arose.
 * I suspect some tests might be necessary in the near or remote future.
 *
 * For now functionality will be covered under STU tests, since there's GPU
 * access there.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE fast_triad_alignment

#include <cstdlib>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "fastsim.h"
#include "fingerprintdb_cuda.h"


using namespace schrodinger::fastsim;
using std::vector;

using schrodinger::fastsim::FastSimServer;

BOOST_AUTO_TEST_CASE(CompareGPUtoCPU)
{
    // TODO:  Replace Schrodinger function that is commented out here
    // Only run this if there's a GPU available
    // if(!schrodinger::gpgpu::is_any_gpu_available()) return;

    FastSimServer server("small.fsim");
    // Fetch a fingerprint to search against, this should always
    // guarantee a 100% match
    const Fingerprint& fp = server.getFingerprint(std::rand() % 20);

    std::vector<char*> gpu_smiles;
    std::vector<char*> gpu_ids;
    std::vector<float> gpu_scores;
    server.similaritySearch(fp, gpu_smiles, gpu_ids, gpu_scores, CalcType::GPU);

    std::vector<char*> cpu_smiles;
    std::vector<char*> cpu_ids;
    std::vector<float> cpu_scores;
    server.similaritySearch(fp, cpu_smiles, cpu_ids, cpu_scores, CalcType::CPU);

    BOOST_CHECK_EQUAL(gpu_smiles.size(), 10);
    for(unsigned int i=0; i<gpu_smiles.size(); i++) {
        BOOST_CHECK_EQUAL(gpu_smiles[i], cpu_smiles[i]);
    }

}

/* 
 * This tests our custom sort function, modified bubble sort to get
 * O(num_required * len(indices))
 */
BOOST_AUTO_TEST_CASE(CPUSort)
{
    vector<int> indices = {0,1,2,3,4,5};
    vector<float> scores = {1,3,2,4,0,7};
    top_results_bubble_sort(indices, scores, 3);

    // Verify that this pushed the top 3 values to the far right of arrays
    BOOST_CHECK_EQUAL(indices[5], 5);
    BOOST_CHECK_EQUAL(scores[5], 7);

    BOOST_CHECK_EQUAL(indices[3], 1);
    BOOST_CHECK_EQUAL(scores[3], 3);

}
