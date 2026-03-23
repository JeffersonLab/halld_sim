#ifndef RANDOMIZED_SDME_H
#define RANDOMIZED_SDME_H

#include <vector>

/**
 * @file randomized_sdme.h
 * @author Hao Li
 * @William & Mary
 * @date 2025-09
 *
 * @brief Generate randomized Spin Density Matrix Elements (SDMEs) from helicity amplitudes.
 *
 * @disclaimer
 * This code is provided "as is", without warranty of any kind, express or implied,
 *
 * -----------------------------------------------------------------------------
 * ## Overview
 * This library provides a function `randomized_sdme()` that generates a random set
 * of 9 spin-density matrix elements (SDMEs) used in vector-meson photoproduction SDME analyses.
 *
 * ## Helicity Matrix Definition
 * We start from the helicity amplitudes
 *
 *   M_{λγ, λV}
 *
 * where
 *   - λγ ∈ {+1, −1} is the photon helicity,
 *   - λV ∈ {+1, 0, −1} is the vector meson helicity.
 *
 * This defines a 2×3 complex matrix of amplitudes:
 *
 *       [ M_{−1,−1}   M_{−1,0}   M_{−1,+1} ]
 *   M = [ M_{+1,−1}   M_{+1,0}   M_{+1,+1} ] .
 *
 * Each matrix element is generated randomly:
 *   - If `use_normal_dist = true`, samples come from a Gaussian(0,1).
 *   - If `use_normal_dist = false`, samples come from a Uniform(−1,1).
 *
 * Both real and imaginary parts of each amplitude are independently randomized.
 *
 * ## Normalization
 * The amplitudes are normalized via
 *
 *   N = (1/2) * Σ |M_{λγ,λV}|²
 *
 * with a small safety offset added to avoid division by zero.
 *
 * ## SDME Construction
 * The spin-density matrix elements are bilinear combinations of the helicity amplitudes,
 * following the definitions in
 *
 *   K. Schilling, P. Seyboth, and G. Wolf, Nucl. Phys. B 15, 397 (1970);
 *   V. Mathieu et al., Phys. Rev. D 97, 094003 (2018).
 *
 * Explicitly, the returned SDMEs are:
 *
 *   ρ_{00}^0     = (1/(2N)) Re[ M_{+1,0}   M_{+1,0}^{*}   + M_{−1,0}   M_{−1,0}^{*}   ]
 *   Re(ρ_{10}^0) = (1/(2N)) Re[ (M_{+1,+1} − M_{+1,−1}) M_{+1,0}^{*} 
 *                             + (M_{−1,+1} − M_{−1,−1}) M_{−1,0}^{*} ]
 *   ρ_{1−1}^0    = (1/(2N)) Re[ M_{+1,+1}  M_{+1,−1}^{*}  + M_{−1,+1}  M_{−1,−1}^{*} ]
 *
 *   ρ_{11}^1     = (1/(2N)) Re[ M_{−1,+1}  M_{+1,+1}^{*}  + M_{+1,−1}  M_{−1,−1}^{*} ]
 *   ρ_{00}^1     = (1/(2N)) Re[ M_{−1,0}   M_{+1,0}^{*}   + M_{+1,0}   M_{−1,0}^{*}  ]
 *
 *   ρ_{1−1}^1    = (1/(2N)) Re[ M_{−1,+1}  M_{+1,−1}^{*}  + M_{+1,+1}  M_{−1,−1}^{*} ]
 *   Im(ρ_{1−1}^2)= (1/(2N)) Im[ M_{−1,+1}  M_{+1,−1}^{*}  − M_{+1,+1}  M_{−1,−1}^{*} ]
 *
 *   Re(ρ_{10}^1) = (1/(2N)) Re[ M_{−1,+1}  M_{+1,0}^{*}   + M_{+1,+1}  M_{−1,0}^{*}  ]
 *   Re(ρ_{10}^2) = (1/(2N)) Re[ M_{−1,+1}  M_{+1,0}^{*}   − M_{+1,+1}  M_{−1,0}^{*}  ]
 *
 * In addition, the sampled SDMEs satisfy the physical constraints, as they are derived
 * in table 2 of Schilling et al. (1970):
  * In addition, the sampled SDMEs satisfy the physical constraints,
 * as derived in Table 2 of Schilling et al. (Nucl. Phys. B15, 397 (1970)):
 *
 *  (Eq. 1)   0 ≤ ρ_{00}^0 ≤ 1
 *
 *  (Eq. 2)   |ρ_{1−1}^0| ≤ (1/2)(1 − ρ_{00}^0)
 *
 *  (Eq. 3)   (Re ρ_{10}^0)^2 ≤ (1/4) ρ_{00}^0 (2 − ρ_{00}^0 − Re ρ_{1−1}^0)
 *
 *  (Eq. 4 not available for Vector meson)   |Im ρ_{1−1}^3| ≤ (1/2)(1 − ρ_{00}^0)
 *
 *  (Eq. 5 not available for Vector meson)   |Im ρ_{10}^3| ≤ sqrt[ (1/2) ρ_{00}^0 (1 − ρ_{00}^0) ]
 *
 *  (Eq. 6)   |ρ_{00}^1| ≤ ρ_{00}^0
 *
 *  (Eq. 7)   |ρ_{11}^1| ≤ (1/2)(1 − ρ_{00}^0)
 *
 *  (Eq. 8)   |ρ_{1−1}^1| ≤ (1/2)(1 − ρ_{00}^0)
 *
 *  (Eq. 9)   |Re ρ_{10}^1| ≤ sqrt[ (1/2) ρ_{00}^0 (1 − ρ_{00}^0) ]
 *
 * (Eq. 10)   |Im ρ_{1−1}^2| ≤ (1/2)(1 − ρ_{00}^0)
 *
 * (Eq. 11)   |Im ρ_{10}^2| ≤ sqrt[ (1/2) ρ_{00}^0 (1 − ρ_{00}^0) ]
 *
 * (Eq. 12)   (Re ρ_{10}^0 ± Re ρ_{10}^1)^2 ≤ (1/8) [ ( (1/2) + (1/2) ρ_{00}^0 ∓ ρ_{11}^1 )
 *                                                   − (ρ_{1−1}^0 ± ρ_{1−1}^1) ]^2
 *                                             − [ ( (3/2) ρ_{00}^0 ± ρ_{00}^1 ∓ ρ_{11}^1
 *                                                   + ρ_{1−1}^0 ± ρ_{1−1}^1 − (1/2) ) ]^2
 *
 * The function returns these 9 real-valued SDMEs in a std::vector<double>.
 *
 * @param use_normal_dist  If true → Gaussian(0,1) distribution.
 *                         If false → Uniform(−1,1) distribution.
 *
 * @return std::vector<double> containing 9 SDMEs.
 */
std::vector<double> randomized_sdme(bool use_normal_dist = false);

#endif // RANDOMIZED_SDME_H
