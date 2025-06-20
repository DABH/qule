#include "GradientDescentPenalty.h"

pair<cx_mat,int> GradientDescentPenalty::Learn(const cx_mat& U_0,const cx_mat& X,const cx_mat& Y) {
    cx_mat U=U_0;
    int counter=0;
	double prevFro = 0;
    while(true||counter<100000) {
        cx_mat dU=objective->DValue(U,X,Y);
        for(unsigned int i=0;i<constraints.size();++i)
            dU+=penalties[i]*constraints[i]->DValue(U,X,Y);
        dU*=learning_rate->Rate();
        U-=dU;

        // if(objective->Value(U,X,Y)<=TOLERANCE)
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

        counter++;
        if(DEBUG_FLAG)
            std::cout<<"counter = "<<counter<<", fro = "<<objective->Value(U,X,Y)<<", mmadU = "<<max(max(abs(dU)))<<std::endl;
    }
    return pair<cx_mat,int>(U,counter);
}
