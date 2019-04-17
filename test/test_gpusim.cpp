#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gpusimilarity

#include <cstdlib>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "calculation_functors.h"
#include "gpusim.h"
#include "fingerprintdb_cuda.h"


using namespace gpusim;
using std::vector;

using gpusim::GPUSimServer;

bool missing_cuda_skip()
{
    static auto cuda_skip = std::getenv("SKIP_CUDA");
    if(cuda_skip != nullptr) return true;
    if(gpusim::get_gpu_count() == 0) return true;

    return false;
}

BOOST_AUTO_TEST_CASE(CompareGPUtoCPU)
{
    if(missing_cuda_skip()) return;

    // TODO:  Replace Schrodinger function that is commented out here
    // Only run this if there's a GPU available
    // if(!schrodinger::gpgpu::is_any_gpu_available()) return;

    QStringList db_fnames;
    db_fnames << "small.fsim";
    GPUSimServer server(db_fnames);
    // Fetch a fingerprint to search against, this should always
    // guarantee a 100% match
    const Fingerprint& fp = server.getFingerprint(std::rand() % 20, "small");
    const vector<int> return_counts = {10, 15};
    const float similarity_cutoff = 0;
    const QString dbkey = "pass"; // Passed to database during creation

    for(auto return_count : return_counts) {
        std::vector<char*> gpu_smiles;
        std::vector<char*> gpu_ids;
        std::vector<float> gpu_scores;
        server.similaritySearch(fp, "small", dbkey, gpu_smiles, gpu_ids,
                gpu_scores, return_count, similarity_cutoff, CalcType::GPU);

        std::vector<char*> cpu_smiles;
        std::vector<char*> cpu_ids;
        std::vector<float> cpu_scores;
        server.similaritySearch(fp, "small", dbkey, cpu_smiles, cpu_ids,
                cpu_scores, return_count, similarity_cutoff, CalcType::CPU);

        BOOST_CHECK_EQUAL(gpu_smiles.size(), return_count);
        for(unsigned int i=0; i<gpu_smiles.size(); i++) {
            BOOST_CHECK_EQUAL(gpu_smiles[i], cpu_smiles[i]);
        }
    }
}


BOOST_AUTO_TEST_CASE(TestSearchMultiple)
{
    if(missing_cuda_skip()) return;
    QStringList db_fnames;
    db_fnames << "small.fsim";
    db_fnames << "small_copy.fsim";
    GPUSimServer server(db_fnames);
    // Fetch a fingerprint to search against, this should always
    // guarantee a 100% match
    const Fingerprint& fp = server.getFingerprint(std::rand() % 20, "small");
    int return_count = 10;
    const float similarity_cutoff = 0;

    std::vector<char*> smiles;
    std::vector<char*> ids;
    std::vector<float> scores;
    std::map<QString, QString> dbname_to_key;
    dbname_to_key["small"] = "pass";
    dbname_to_key["small_copy"] = "pass";

    server.searchDatabases(fp, return_count, similarity_cutoff, dbname_to_key,
            smiles, ids, scores);

    BOOST_REQUIRE_EQUAL(smiles.size(), return_count);
    // Two copies of the database should always have duplicate top 2 results
    BOOST_REQUIRE_EQUAL(smiles[0], smiles[1]);
}


BOOST_AUTO_TEST_CASE(TestSimilarityCutoff)
{
    if(missing_cuda_skip()) return;
    QStringList db_fnames;
    db_fnames << "small.fsim";
    GPUSimServer server(db_fnames);
    const Fingerprint& fp = server.getFingerprint(0, "small");
    int return_count = 10;
    const vector<float> similarity_cutoffs = {0, 0.1, 0.3, 0.4};
    // Results returned for each cutoff range
    const vector<int> result_counts = {10, 10, 3, 1};
    const QString dbkey = "pass"; // Passed to database during creation

    for(unsigned int i=0; i<similarity_cutoffs.size(); i++) {
        std::vector<char*> smiles;
        std::vector<char*> ids;
        std::vector<float> scores;
        server.similaritySearch(fp, "small", dbkey, smiles, ids, scores,
                return_count, similarity_cutoffs[i], CalcType::GPU);

        BOOST_CHECK_EQUAL(smiles.size(), result_counts[i]);
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
    BOOST_CHECK_EQUAL(indices[0], 5);
    BOOST_CHECK_EQUAL(scores[0], 7);

    BOOST_CHECK_EQUAL(indices[2], 1);
    BOOST_CHECK_EQUAL(scores[2], 3);

}

BOOST_AUTO_TEST_CASE(FoldFingerprint)
{
    int factor = 2;
    vector<int> fp = {32, 24, 11, 7};
    vector<int> ref_answer = {43, 31};
    vector<int> answer(fp.size() / factor);
    gpusim::FoldFingerprintFunctorCPU(factor, fp.size(), fp, answer)(0);

    for(unsigned int i=0; i<ref_answer.size(); i++){
    	BOOST_CHECK_EQUAL(answer[i], ref_answer[i]);
    }

    factor = 4;
    answer.resize(1);
    answer[0] = 0;
    gpusim::FoldFingerprintFunctorCPU(factor, fp.size(), fp, answer)(0);
    BOOST_CHECK_EQUAL(answer.size(), 1);
    BOOST_CHECK_EQUAL(answer[0], 63);
}

BOOST_AUTO_TEST_CASE(getNextGPU)
{
    if(missing_cuda_skip()) return;
    unsigned int gpucount = gpusim::get_gpu_count();
    // Make sure first go around gets every valid GPU
    for(unsigned int i=0; i<gpucount; i++) {
        BOOST_CHECK_EQUAL(i, gpusim::get_next_gpu(1));
    }
    // Make sure second go around loops through them all again
    for(unsigned int i=0; i<gpucount; i++) {
        BOOST_CHECK_EQUAL(i, gpusim::get_next_gpu(1));
    }
}
