#include "DeltaLearningPenalty.h"

pair<cx_mat, int> DeltaLearningPenalty::Learn(const cx_mat &U_0, const cx_mat &X, const cx_mat &Y)
{
	cx_mat U = U_0;
	int counter = 0;
	double prevFro = 0;

	while (true || counter < 10)
	{
		cx_mat dU;
		for (unsigned int m = 0; m < X.n_cols; ++m)
		{
			dU = objective->DValue(U, X.col(m), Y.col(m));
			for (unsigned int i = 0; i < constraints.size(); ++i)
				dU += penalties[i] * constraints[i]->DValue(U, X.col(m), Y.col(m));
			dU *= learning_rate->Rate();
			U -= dU;
		}

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

		// if (objective->Value(U, X, Y) < TOLERANCE) //{
		// 	// std::cout<<"LAST OBJECTIVE VALUE IS "<<objective->Value(U,X,Y)<<", TOLERANCE is "<<TOLERANCE<<std::endl;
		// 	// break;
		// 	//}
		// 	// if(max(max(abs(dU)))<TOLERANCE||objective->Value(U,X,Y)<TOLERANCE)
		// 	// if(objective->Value(U,X,Y)<=1.e-16)//max(max(abs(dU)))<1e-16)//objective->Value(U,X,Y)<=1.e-16)
		// 	break;
		// // std::cout<<U<<std::endl;
		counter++;
		if (DEBUG_FLAG)
			std::cout << "counter = " << counter << ", fro = " << objective->Value(U, X, Y) << ", mmadU = " << max(max(abs(dU))) << std::endl;
	}
	return pair<cx_mat, int>(U, counter);
}
