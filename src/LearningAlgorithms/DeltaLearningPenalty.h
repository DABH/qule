#ifndef DELTA_LEARNING_PENALTY_H
#define DELTA_LEARNING_PENALTY_H

#include "LearningAlgorithm.h"
#include "../Globals.h"

using namespace arma;
using std::vector;

class DeltaLearningPenalty : public LearningAlgorithm {
    public:
        DeltaLearningPenalty(const Function* _objective,const vector<Function*>& _constraints,const LearningRate* _learning_rate,const vector<double>& _penalties): LearningAlgorithm(_objective,_constraints,_learning_rate),penalties(_penalties){} //TODO: Add a check to ensure that the number of penalties matches the number of constraints
        pair<cx_mat,int> Learn(const cx_mat& U_0,const cx_mat& X,const cx_mat& Y);
    private:
        vector<double> penalties;
};

#endif
