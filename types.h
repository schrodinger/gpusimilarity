#ifndef FP_TYPES_H
#define FP_TYPES_H

#include <vector>
#include <utility>

namespace fastsim
{

typedef std::pair<std::vector<char*>, std::vector<float> > SimResults;
typedef std::vector<int> Fingerprint;
}

#endif
