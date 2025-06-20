#include "Orthogonalization.h"

double Orthogonalization::Value(const cx_mat& U,const cx_mat&,const cx_mat&) const {
    return 0.25*pow(norm(U.t()*U-eye<cx_mat>(U.n_cols,U.n_cols),"fro"),2.0);
}

cx_mat Orthogonalization::DValue(const cx_mat& U,const cx_mat&,const cx_mat&) const {
    return U*(U.t()*U-eye<cx_mat>(U.n_cols,U.n_cols));
}

cx_mat Orthogonalization::DDValue(const cx_mat&,const cx_mat&,const cx_mat&) const{
    return cx_mat();
}
