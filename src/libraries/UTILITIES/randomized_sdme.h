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
 *   V. Mathieu et al., Phys. Rev. D 97, 094003 (2018).
 *
 * Explicitly, the returned SDMEs are:
 *
 *     ρ_{00}^0   = (1/N) Re[ M_{+1,0} M_{+1,0}* ]
 *  Re(ρ_{10}^0)  = (1/(2N)) Re[ (M_{+1,+1} − M_{+1,−1}) M_{+1,0}* ] 
 *     ρ_{1−1}^0  = (1/N) Re[ M_{+1,+1} M_{+1,−1}* ]
 *     ρ_{11}^1   = (1/N) Re[ M_{−1,+1} M_{+1,+1}* ]
 *     ρ_{00}^1   = (1/N) Re[ M_{−1,0}  M_{+1,0}* ]
 *     ρ_{1−1}^1  = (1/N) Re[ M_{−1,+1} M_{+1,−1}* + M_{+1,+1} M_{−1,−1}* ]
 *  Im(ρ_{1−1}^2) = (1/N) Im[ M_{−1,+1} M_{+1,−1}* − M_{+1,+1} M_{−1,−1}* ] 
 *  Re(ρ_{10}^1)  = (1/N) Re[ M_{−1,+1} M_{+1,0}*  + M_{+1,+1} M_{−1,0}*  ]
 *  Re(ρ_{10}^2)  = (1/N) Re[ M_{−1,+1} M_{+1,0}*  − M_{+1,+1} M_{−1,0}*  ]

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
