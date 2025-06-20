#ifndef CONSTANT_H
#define CONSTANT_H

#include "LearningRate.h"

using namespace arma;
using std::vector;

class Constant : public LearningRate {
    public:
        Constant(const double & _rate): rate(_rate){};
        double Rate() const;
        double Rate(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk) const;
    private:
        double rate;
};

#endif
