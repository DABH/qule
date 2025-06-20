#ifndef FUNCTION_H
#define FUNCTION_H

#include <armadillo>

using namespace arma;

class Function {
public:
    Function(){};

    virtual double Value(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const=0;

    // returns a cx_mat of the same size as U
    virtual cx_mat DValue(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const=0;

    // returns a cx_mat of the size n^2 * n^2 (where n is size of U)
    virtual cx_mat DDValue(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const=0;
};

#endif
