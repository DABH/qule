#ifndef NEWTONS_METHOD_H
#define NEWTONS_METHOD_H

// \warning: Newton's Method is designed solely for unconstrained Procrustes Problem. Will need to work out Hessian of Unitarization constraint

#include "LearningAlgorithm.h"

using namespace arma;
using std::vector;

class NewtonsMethod : public LearningAlgorithm {
    public:
        NewtonsMethod(const Function* _objective,const std::vector<Function*>& _constraints,const LearningRate* _learning_rate): LearningAlgorithm(_objective,_constraints,_learning_rate){}
        pair<cx_mat,int> Learn(const cx_mat& U_0,const cx_mat& X,const cx_mat& Y);
    private:
        cx_vec MatToVec(const cx_mat& X) const;
};

#endif
