#ifndef FP_TYPES_H
#define FP_TYPES_H

#include <vector>
#include <utility>

namespace fastsim
{

const int RETURN_COUNT = 10;

typedef std::pair<std::vector<char*>, std::vector<float> > SimResults;
typedef std::vector<int> Fingerprint;
}

#endif
