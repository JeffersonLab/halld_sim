#include "randomized_sdme.h"
#include <random>
#include <complex>

/**
 * Hao Li
 * @William & Mary
 * 2025-09
 *
 * Generate randomized Spin Density Matrix Elements (SDMEs) from helicity amplitudes.
 * -----------------------------------------------------------------------------
 */

std::vector<double> randomized_sdme(bool use_normal_dist) {
    // Random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
  
    // Distributions
    static std::normal_distribution<double> normal_dist(0.0, 1.0);
    static std::uniform_real_distribution<double> uniform_dist(-1.0, 1.0);
  
    // Unified sampling function
    auto sample = [&]() {
        return use_normal_dist ? normal_dist(gen) : uniform_dist(gen);
    };
  
    // generate 2x3 complex matrix:
    //    helicity amplitudes M_{lambda_g, lambda_V},
    //    where lambda_g = 1, -1 and lambda_V = 1, 0, -1,
    //    then compute normalization N = 0.5 * sum |M|^2.
    std::complex<double> M[2][3];
    for (int i = 0; i < 2; ++i){
      for (int j = 0; j < 3; ++j){
        M[i][j] = std::complex<double>(sample(), sample());
      }
    }
        
  
    // human readable shortcuts 
    const std::complex<double>& M_m1_m1 = M[0][0];  // M_{-1, -1}
    const std::complex<double>& M_p1_m1 = M[1][0];  // M_{+1, -1}
    const std::complex<double>& M_m1_0  = M[0][1];  // M_{-1, 0}
    const std::complex<double>& M_p1_0  = M[1][1];  // M_{+1, 0}
    const std::complex<double>& M_m1_p1 = M[0][2];  // M_{-1, +1}
    const std::complex<double>& M_p1_p1 = M[1][2];  // M_{+1, +1}
  
    // compute SDMEs, from [Schilling; MATHIEU et al. PHYS. REV. D 97, 094003]
    const double N =
    0.5 * ( std::norm(M_p1_m1) + std::norm(M_p1_0) + std::norm(M_p1_p1)
          + std::norm(M_m1_m1) + std::norm(M_m1_0) + std::norm(M_m1_p1) );
    double Ninv = 1.0 / (1e-12 + N); // small offset to avoid division by zero
    
    // (B1a-c)
    double rho000  = 0.5 * Ninv * ( (M_p1_0 * std::conj(M_p1_0)).real()
                                  + (M_m1_0 * std::conj(M_m1_0)).real() );

    double rho100  = 0.5 * Ninv * ( ((M_p1_p1 - M_p1_m1) * std::conj(M_p1_0)).real()
                                  + ((M_m1_p1 - M_m1_m1) * std::conj(M_m1_0)).real() );

    double rho1m10 = 0.5 * Ninv * ( (M_p1_p1 * std::conj(M_p1_m1)).real()
                                  + (M_m1_p1 * std::conj(M_m1_m1)).real() );         
    
    // (B1d-e)
    double rho111  = 0.5 * Ninv * (M_m1_p1 * std::conj(M_p1_p1) + M_p1_m1 * std::conj(M_m1_m1)).real();
    double rho001  = 0.5 * Ninv * (M_m1_0  * std::conj(M_p1_0 ) + M_p1_0  * std::conj(M_m1_0 )).real();
  
    // (B1f-g)
    std::complex<double> prod_m1p1 = M_m1_p1 * std::conj(M_p1_m1);
    std::complex<double> prod_11   = M_p1_p1 * std::conj(M_m1_m1);
    double rho1m11 = 0.5 * Ninv * (prod_m1p1 + prod_11).real();
    double rho1m12 = 0.5 * Ninv * (prod_m1p1 - prod_11).imag();
  
    // (B1h-i)
    double re_prod_10   = (M_m1_p1 * std::conj(M_p1_0)).real();
    double re_prod_m10  = (M_p1_p1 * std::conj(M_m1_0)).real();
    double rho101 = 0.5 * Ninv * (re_prod_10 + re_prod_m10);
    double rho102 = 0.5 * Ninv * (re_prod_10 - re_prod_m10);
  
    // Output vector of 9 SDMEs
    return {
        rho000, rho100, rho1m10,
        rho111, rho001, rho101,
        rho1m11, rho102, rho1m12
    };
}
  
