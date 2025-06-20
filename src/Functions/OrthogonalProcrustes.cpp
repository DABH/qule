#include "OrthogonalProcrustes.h"

double OrthogonalProcrustes::Value(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const {
    return 0.5*pow(norm(U*X-Y,"fro"),2.0);
}

cx_mat OrthogonalProcrustes::DValue(const cx_mat& U,const cx_mat& X,const cx_mat& Y) const {
    return (U*X-Y)*X.t();
}

cx_mat OrthogonalProcrustes::DDValue(const cx_mat& U,const cx_mat& X,const cx_mat&) const {
    cx_mat XHX=X*X.t(); //Allocate Hermitian of X beforehand
    // Using indexing H_{i,j,k,l} = H_{i*n_cols+j,k*n_cols+l}
    cx_mat H(U.n_rows*U.n_cols,U.n_rows*U.n_cols,fill::zeros);
    for(unsigned int i=0;i<U.n_rows;++i)
        for(unsigned int j=0;j<U.n_cols;++j)
            for(unsigned int k=0;k<U.n_rows;++k)
                for(unsigned int l=0;l<U.n_cols;++l)
                    if(i==k)
                        H(i*U.n_cols+j,k*U.n_cols+l)=XHX(l,j);
    return H;
}
