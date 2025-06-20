#include "NewtonsMethod.h"
#include "../Globals.h"

// Turn X_{ij} into X_{i*numcols+j} - CAN PARALLELIZE THIS
cx_vec NewtonsMethod::MatToVec(const cx_mat& X) const {
    cx_vec Y(X.n_cols*X.n_rows);
    for(unsigned int i=0;i<X.n_rows;++i)
        for(unsigned int j=0;j<X.n_cols;++j)
            Y(i*X.n_cols+j)=X(i,j);
    return Y;
}

pair<cx_mat,int> NewtonsMethod::Learn(const cx_mat& U_0,const cx_mat& X,const cx_mat& Y) {
    cx_mat U=U_0;
    int counter=0;
    double prevFro = 0;
    while(true||counter<10) {
        cx_mat dU=objective->DValue(U,X,Y);
        cx_vec dUvec=MatToVec(dU);
        cx_mat H=objective->DDValue(U,X,Y);
        cx_vec Solvec;
        solve(Solvec,H,dUvec);
        double alpha=learning_rate->Rate();
        //double alpha=learning_rate->Rate(U,X,Y,dU);
        //std::cout<<"alpha = "<<alpha<<std::endl;
        Solvec*=alpha;
        for(unsigned int i=0;i<U.n_rows;++i)
            for(unsigned int j=0;j<U.n_cols;++j)
                U(i,j)-=Solvec(i*U.n_cols+j); // PARALLELIZE THIS
        counter++;
        if(DEBUG_FLAG)
            std::cout<<"counter = "<<counter<<", obj = "<<objective->Value(U,X,Y)<<", mmadU = "<<max(max(abs(dU)))<<std::endl;
        //if(objective->Value(U,X,Y)<=TOLERANCE)
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
        //std::cout<<U<<std::endl;
    }
    return pair<cx_mat,int>(U,counter);
}
