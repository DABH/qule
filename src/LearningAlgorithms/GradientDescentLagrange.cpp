#include "GradientDescentLagrange.h"

pair<cx_mat, int> GradientDescentLagrange::Learn(const cx_mat &U_0, const cx_mat &X, const cx_mat &Y)
{
	cx_mat U = U_0;
	std::vector<double> lambda(constraints.size(), 1.0);
	int counter = 0;
	double prevFro = 0;

	// why is this while loop true || counter < 100?
	while (true || counter < 100)
	{
		// std::cout<<U<<std::endl;
		counter++;
		if (DEBUG_FLAG)
			std::cout << "counter = " << counter << ", fro = " << objective->Value(U, X, Y) << std::endl;

		cx_mat dU = objective->DValue(U, X, Y);
		for (unsigned int i = 0; i < constraints.size(); ++i)
		{
			// below is where we adjust whether we make it it + lambda or - lambda
			// since U -= dU, if we make du += lambda that is the same as - lambda, and -= lambda -> +
			// std::cout << lambda[i]*constraints[i]->DValue(U,X,Y) << std::endl;

			dU += lambda[i] * constraints[i]->DValue(U, X, Y);
			lambda[i] += learning_rate->Rate() * constraints[i]->Value(U, X, Y);
		}

		dU *= learning_rate->Rate();
		U -= dU;

        int check = checkTolerance(objective->Value(U,X,Y), prevFro);
        if(check == 1) // converged
			break;
        else if (check == -1) // stalled
        {
            counter = -1;
            break;
        }
		else // neither
			prevFro = objective->Value(U,X,Y);
	}
	return pair<cx_mat, int>(U, counter);
}
