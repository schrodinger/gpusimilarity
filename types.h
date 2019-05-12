#ifndef FP_TYPES_H
#define FP_TYPES_H

#include <utility>
#include <vector>

namespace gpusim
{

typedef std::pair<std::vector<char*>, std::vector<float>> SimResults;
typedef std::vector<int> Fingerprint;
} // namespace gpusim

#endif
