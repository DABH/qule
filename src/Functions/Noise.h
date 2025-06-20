#ifndef NOISE_H
#define NOISE_H

#include "Function.h"

using namespace arma;

class GaussianNoise {

public:
    GaussianNoise(){};

	cx_mat addNoise(const cx_mat& X, double magnitude) const;
};

#endif
