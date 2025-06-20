#ifndef ORTHOGONAL_PROCRUSTES_H
#define ORTHOGONAL_PROCRUSTES_H

#include "Function.h"

using namespace arma;

class OrthogonalProcrustes : public Function {

public:
    OrthogonalProcrustes(){};

    double Value(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const;

    // returns a cx_mat of the same size as U
    cx_mat DValue(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const;

    // returns a cx_mat of the size n^2 * n^2 (where n is size of U)
    cx_mat DDValue(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const;
};

#endif
