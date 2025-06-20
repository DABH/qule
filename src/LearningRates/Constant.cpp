#include "Constant.h"

double Constant::Rate() const {
    return rate;
}

double Constant::Rate(const cx_mat&,const cx_mat&,const cx_mat&,const cx_mat&) const {
    return rate;
}
