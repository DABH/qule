#ifndef ADAM_H
#define ADAM_H

#include "LearningRate.h"
#include <armadillo>
#include <fstream>

class Adam : public LearningRate {
public:
    // The "log_stream" parameter is first and mandatory,
    // but the other four have defaults.
    Adam(
         std::ofstream &log_stream,
         double learning_rate = 0.1, 
         double beta1         = 0.9, 
         double beta2         = 0.999, 
         double epsilon       = 1e-8   // <-- NO trailing comma
    );

    virtual double Rate(const arma::cx_mat& Uk, 
                        const arma::cx_mat& X, 
                        const arma::cx_mat& Y, 
                        const arma::cx_mat& dUk) const override;

private:
    // Put these in the same order you'll initialize them in Adam.cpp:
    // 1) log_stream_   (since it appears first in the constructor's initializer list)
    // 2) learning_rate
    // 3) beta1, beta2, epsilon
    // 4) t
    // 5) m, v
    // This ensures no reorder warnings.

    std::ofstream &log_stream_;

    double learning_rate;  
    double beta1;
    double beta2;
    double epsilon;

    mutable int t;         // Timestep counter
    mutable arma::cx_mat m, v;  // Moment vectors
};

#endif