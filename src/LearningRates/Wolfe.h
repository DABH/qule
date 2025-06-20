#ifndef WOLFE_H
#define WOLFE_H

#include "LearningRate.h"
#include <fstream>  

using namespace arma;
using std::vector;

class Wolfe : public LearningRate {
    public:
        // Wolfe(){}
        // Let the constructor take a reference (or pointer) to the log stream
        Wolfe(std::ofstream &logStream) 
        : logStream_(logStream) {}
        double Rate(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk) const;
    //private:
        //double rate;
    private:
        // reference to the output stream
        std::ofstream &logStream_;
        double phi(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha) const;
        double phi_prime(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha) const;
        double zoom(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha_left,double alpha_right) const;
        double c1=1e-4; // values from Nocedal and Wright p. 62
        double c2=0.9;
};

#endif
