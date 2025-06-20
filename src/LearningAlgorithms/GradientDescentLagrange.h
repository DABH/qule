#ifndef GRADIENT_DESCENT_LAGRANGE_H
#define GRADIENT_DESCENT_LAGRANGE_H

#include "LearningAlgorithm.h"

using namespace arma;
using std::vector;

class GradientDescentLagrange : public LearningAlgorithm {
    public:
        GradientDescentLagrange(const Function* _objective,const vector<Function*>& _constraints,const LearningRate* _learning_rate): LearningAlgorithm(_objective,_constraints,_learning_rate){}
        pair<cx_mat,int> Learn(const cx_mat& U_0,const cx_mat& X,const cx_mat& Y);
};

#endif
