#ifndef LEARNING_RATE_H
#define LEARNING_RATE_H

#include <armadillo>
#include <vector>

using namespace arma;
using std::vector;

class LearningRate {
    public:
        LearningRate(){};
        virtual double Rate() const { return 0.0; }; // fill this in later
        virtual double Rate(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk) const = 0;// { return 0.0; }; // fill this in later
};

#endif
