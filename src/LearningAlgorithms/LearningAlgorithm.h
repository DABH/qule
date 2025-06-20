#ifndef LEARNING_ALGORITHM_H
#define LEARNING_ALGORITHM_H

#include <armadillo>
#include <vector>
#include "../Globals.h"
#include "../Functions/Function.h"
#include "../LearningRates/LearningRate.h"

#define LEARNING_ALGO_NUM_OF_ZERO_TO_EXIT 30000

using namespace arma;
using std::pair;
using std::vector;

class LearningAlgorithm
{
public:
	LearningAlgorithm(const Function *_objective, const vector<Function *> &_constraints, const LearningRate *_learning_rate) : objective(_objective), constraints(_constraints), learning_rate(_learning_rate) {}
	virtual pair<cx_mat, int> Learn(const cx_mat &U_0, const cx_mat &X, const cx_mat &Y) = 0;
	int checkTolerance(const double currVal, const double prevVal)
	{
		// std::cout << "Checking Tolerance: Curr Val and Prev Val" << currVal << " " << prevVal << std::endl;
		// std::cout << (double)(currVal - prevVal) << std::endl;
		double diffInVal = (double)(currVal - prevVal);
		if (diffInVal < 0)
			diffInVal = -1 * diffInVal;

		// sometimes the difference between two steps is 0 (floating point precision?)
		// but we don't want to just yet exit out, only want to exit out if it has been zero for a long time.
		// if (diffInVal == 0)
		// {
		// 	std::cout << "0 diff" << std::endl;
		// 	countOfZerosInTolerance++;
		// 	if (this->countOfZerosInTolerance > LEARNING_ALGO_NUM_OF_ZERO_TO_EXIT)
		// 	{
		// 		return true;
		// 	}
		// 	return currVal < TOLERANCE;
		// }
		// else
		// {
		if (currVal < TOLERANCE)
			return 1;
		else if (diffInVal < FRO_CHANGE_TOLERANCE)
			return -1;
		else
			return 0;
		// }
	};

protected:
	const Function *objective;
	const vector<Function *> constraints;
	const LearningRate *learning_rate;
	int countOfZerosInTolerance = 0;
};

#endif
