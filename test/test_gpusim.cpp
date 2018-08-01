/*
 * NOTE:  I mainly left this in here because this is working with waf, and
 * I didn't want to have to figure out how to do this again if problems arose.
 * I suspect some tests might be necessary in the near or remote future.
 *
 * For now functionality will be covered under STU tests, since there's GPU
 * access there.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gpusimilarity

#include <cstdlib>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "gpusim.h"
#include "fingerprintdb_cuda.h"


using namespace gpusim;
using std::vector;

using gpusim::GPUSimServer;

BOOST_AUTO_TEST_CASE(CompareGPUtoCPU)
{
    // TODO:  Replace Schrodinger function that is commented out here
    // Only run this if there's a GPU available
    // if(!schrodinger::gpgpu::is_any_gpu_available()) return;

    QStringList db_fnames;
    db_fnames << "small.fsim";
    GPUSimServer server(db_fnames);
    // Fetch a fingerprint to search against, this should always
    // guarantee a 100% match
    const Fingerprint& fp = server.getFingerprint(std::rand() % 20);
    const vector<int> return_counts = {10, 15};

    for(auto return_count : return_counts) {
        std::vector<char*> gpu_smiles;
        std::vector<char*> gpu_ids;
        std::vector<float> gpu_scores;
        server.similaritySearch(fp, gpu_smiles, gpu_ids, gpu_scores, return_count, CalcType::GPU);

        std::vector<char*> cpu_smiles;
        std::vector<char*> cpu_ids;
        std::vector<float> cpu_scores;
        server.similaritySearch(fp, cpu_smiles, cpu_ids, cpu_scores, return_count, CalcType::CPU);

        BOOST_CHECK_EQUAL(gpu_smiles.size(), return_count);
        for(unsigned int i=0; i<gpu_smiles.size(); i++) {
            BOOST_CHECK_EQUAL(gpu_smiles[i], cpu_smiles[i]);
        }
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

BOOST_AUTO_TEST_CASE(FoldFingerprint)
{
    int factor = 2;
    vector<int> fp = {32, 24, 11, 7};
    vector<int> ref_answer = {43, 31};
    vector<int> answer = fold_fingerprint(fp, factor);
    for(int i=0; i<ref_answer.size(); i++){
    	BOOST_CHECK_EQUAL(answer[i], ref_answer[i]);
    }

    factor = 4;
    answer = fold_fingerprint(fp, factor);
    BOOST_CHECK_EQUAL(answer.size(), 1);
    BOOST_CHECK_EQUAL(answer[0], 63);
}
