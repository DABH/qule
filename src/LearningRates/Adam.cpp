#include "Adam.h"

Adam::Adam(std::ofstream &log_stream,
           double learning_rate,
           double beta1,
           double beta2,
           double epsilon
)
// Adam::Adam(
//     double learning_rate,
//     double beta1,
//     double beta2,
//     double epsilon
// )
  // Match the order we declared them in the header:
  : log_stream_(log_stream)
  , learning_rate(learning_rate)
  , beta1(beta1)
  , beta2(beta2)
  , epsilon(epsilon)
  , t(0)
  , m()
  , v()
{
    // Initialize moment vectors
    m.zeros();
    v.zeros();
}

double Adam::Rate(const arma::cx_mat& Uk, 
                  const arma::cx_mat& X, 
                  const arma::cx_mat& Y, 
                  const arma::cx_mat& dUk) const
{
    if (t == 0) {
        m.zeros(dUk.n_rows, dUk.n_cols);
        v.zeros(dUk.n_rows, dUk.n_cols);
    }
    t++;

    // Update biased moment estimates
    m = beta1 * m + (1.0 - beta1) * dUk;
    v = beta2 * v + (1.0 - beta2) * arma::square(dUk);

    // Bias correction
    arma::cx_mat m_hat = m / (1.0 - std::pow(beta1, t));
    arma::cx_mat v_hat = v / (1.0 - std::pow(beta2, t));

    // Compute the step size (real part)
    // double alpha = std::real(
    //     arma::as_scalar(learning_rate * m_hat / (arma::sqrt(v_hat) + epsilon))
    // );
    // Take average across all matrix elements
    // arma::cx_mat stepMat = learning_rate * m_hat / (arma::sqrt(v_hat) + epsilon);
    // auto sumVal = arma::accu(stepMat);            // sum all complex entries
    // double alpha = std::real(sumVal) / stepMat.n_elem;

    // maximum absolute
    arma::cx_mat stepMat = learning_rate * m_hat / (arma::sqrt(v_hat) + epsilon);
    double alpha = stepMat.empty()
        ? 0.0
        : arma::max( arma::max( arma::abs(stepMat) ) );

    
    alpha = std::max(alpha, 0.1);
    alpha = std::min(alpha, 0.35);
    // Log the alpha
    log_stream_ << alpha << "\n";

    return alpha;
}
