#ifndef FP_TYPES_H
#define FP_TYPES_H

#include <vector>
#include <utility>

namespace fastsim
{

const int RETURN_COUNT = 10;

using SimResults = std::pair<std::vector<char*>, std::vector<float>>;
using Fingerprint = std::vector<int>;
}

#endif
