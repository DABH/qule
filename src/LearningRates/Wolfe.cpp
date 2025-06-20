#include "Wolfe.h"
#include <iostream>
// double Wolfe::Rate(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk) const {
//     //return 0.01;
//     double alpha_im1=0;
//     double alpha_max=100000.0;
//     double alpha_i=0.1*alpha_max;
//     int i=1;
//     while(true){
//         if (phi(Uk,X,Y,dUk,alpha_i) > phi(Uk,X,Y,dUk,0) + c1 * alpha_i * phi_prime(Uk,X,Y,dUk,0) ||(phi(Uk,X,Y,dUk,alpha_i) >= phi(Uk,X,Y,dUk,alpha_im1) && i > 1))
//             return zoom(Uk,X,Y,dUk,alpha_im1,alpha_i);
//         if(abs(phi_prime(Uk,X,Y,dUk,alpha_i))<=-c2*phi_prime(Uk,X,Y,dUk,0))
//             return alpha_i;
//         if(phi_prime(Uk,X,Y,dUk,alpha_i)>=0)
//             return zoom(Uk,X,Y,dUk,alpha_i,alpha_im1);
//         alpha_im1=alpha_i;
//         alpha_i=0.5*(alpha_max+alpha_i);
//         std::cout<<"alpha_i = "<<alpha_i<<", phi_prime = "<<phi_prime(Uk,X,Y,dUk,alpha_i)<<", blah = "<<-c2*phi_prime(Uk,X,Y,dUk,0)<<std::endl;
//         i++;
//     }
// }

double Wolfe::Rate(const cx_mat& Uk, const cx_mat& X, const cx_mat& Y, const cx_mat& dUk) const {
    double alpha_im1 = 0;
    double alpha_max = 100000.0;
    double alpha_i = std::min(1.0, 0.1 * alpha_max / norm(dUk, "fro"));
    int i = 1;
    
    while (i < 1000) { // Add iteration limit
        if (phi(Uk, X, Y, dUk, alpha_i) > phi(Uk, X, Y, dUk, 0) + c1 * alpha_i * phi_prime(Uk, X, Y, dUk, 0) ||
            (phi(Uk, X, Y, dUk, alpha_i) >= phi(Uk, X, Y, dUk, alpha_im1) && i > 1)) {
                logStream_ << alpha_i << "\n";
            return zoom(Uk, X, Y, dUk, alpha_im1, alpha_i);
        }
        if (abs(phi_prime(Uk, X, Y, dUk, alpha_i)) <= -c2 * phi_prime(Uk, X, Y, dUk, 0)) {
            logStream_ << alpha_i << "\n";
            return alpha_i;
        }
        if (phi_prime(Uk, X, Y, dUk, alpha_i) >= 0) {
            logStream_ << alpha_i << "\n";
            return zoom(Uk, X, Y, dUk, alpha_i, alpha_im1);
        }
        alpha_im1 = alpha_i;
        alpha_i = std::max(0.5 * (alpha_max + alpha_i), 1e-8);
        i++;
    }
    logStream_ << alpha_i << "\n";
    return alpha_i;
}

double Wolfe::phi(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha) const {
    return 0.5*pow(norm((Uk-alpha*dUk)*X-Y,"fro"),2.0);
}

double Wolfe::phi_prime(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha) const {
    return alpha*real(trace(dUk*X*X.t()*dUk.t()))+trace(real(Y*X.t()*dUk.t()))-trace(real(dUk*X*X.t()*Uk.t()));
    //return (alpha-1.0)*pow(norm(dUk*X,"fro"),2.0)+trace(real(Y*X.t()*dUk.t()));
}

// double Wolfe::zoom(const cx_mat& Uk,const cx_mat& X,const cx_mat& Y,const cx_mat& dUk,double alpha_left,double alpha_right) const {
//     std::cout<<"inside zoom!"<<std::endl;
//     double alpha_lo=alpha_left;
//     double alpha_hi=alpha_right;
//     while(true){
//         std::cout<<"alpha_lo = "<<alpha_lo<<", alpha_hi = "<<alpha_hi<<std::endl;
//         double alpha_j=0.5*(alpha_hi+alpha_lo);
//         //if(fabs(alpha_j-alpha_hi)<=1e-15)return alpha_j; // hack
//         if(phi(Uk,X,Y,dUk,alpha_j)>phi(Uk,X,Y,dUk,0)+c1*alpha_j*phi_prime(Uk,X,Y,dUk,0)||phi(Uk,X,Y,dUk,alpha_j)>=phi(Uk,X,Y,dUk,alpha_lo))
//             alpha_hi=alpha_j;
//         else{
//             if(abs(phi_prime(Uk,X,Y,dUk,alpha_j))<=-c2*phi_prime(Uk,X,Y,dUk,0))
//                 return alpha_j;
//             else
//                 std::cout<<"one = "<<abs(phi_prime(Uk,X,Y,dUk,alpha_j))<<", two = "<<(-c2*phi_prime(Uk,X,Y,dUk,0))<<std::endl;
//             if(phi_prime(Uk,X,Y,dUk,alpha_j)*(alpha_hi-alpha_lo)>=0)
//                 alpha_hi=alpha_lo;
//             alpha_lo=alpha_j;
//         }
//     }
// }

double Wolfe::zoom(const cx_mat& Uk, const cx_mat& X, const cx_mat& Y, const cx_mat& dUk, double alpha_left, double alpha_right) const {
    // std::cout << "inside zoom!" << std::endl;
    double alpha_lo = alpha_left;
    double alpha_hi = alpha_right;
    int zoom_iter = 0;
    int max_zoom_iters = 1000;
    
    while (zoom_iter++ < max_zoom_iters) {
        double alpha_j = 0.5 * (alpha_hi + alpha_lo);
        
        if (std::abs(alpha_hi - alpha_lo) < 1e-8) {
            return alpha_j;
        }
        
        if (phi(Uk, X, Y, dUk, alpha_j) > phi(Uk, X, Y, dUk, 0) + c1 * alpha_j * phi_prime(Uk, X, Y, dUk, 0) ||
            phi(Uk, X, Y, dUk, alpha_j) >= phi(Uk, X, Y, dUk, alpha_lo)) {
            alpha_hi = alpha_j;
        } else {
            if (abs(phi_prime(Uk, X, Y, dUk, alpha_j)) <= -c2 * phi_prime(Uk, X, Y, dUk, 0)) {
                return alpha_j;
            }
            if (phi_prime(Uk, X, Y, dUk, alpha_j) * (alpha_hi - alpha_lo) >= 0) {
                alpha_hi = alpha_lo;
            }
            alpha_lo = alpha_j;
        }
    }

    return (alpha_hi + alpha_lo) / 2;  // Return midpoint as fallback
}