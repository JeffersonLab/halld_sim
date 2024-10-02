#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <math.h>

#include "TLorentzVector.h"

#include "IUAmpTools/Kinematics.h"

#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF0.h"

/* 
 *  Usage: KopfKMatrixF0 <daughter 1> <daughter 2> <channel> <Re[f0(500)]> <Im[f0(500)]> ...
 * 
 *  Parameters:
 *  1. <daughter 1>: a single digit number or pair of digits corresponding to the location of
 *     the first decay product in the kinematic arrays. Note that this array typically starts
 *     with the beam in location 0, so these values could be considered to be 1-indexed. If
 *     the number is more than one digit long, each single digit will represent a particle
 *     which, when added together, form the first decay product (e.g. 23 would sum particles
 *     2 and 3 to make a new state).
 *  2. <daughter 2>: Same as (1) but for the second decay product.
 *  3. <channel>: A single digit corresponding to the final state channel:
 *     0: pi-pi
 *     1: 2pi-2pi (effective channel for all decays into four pions not present in other channels)
 *     2: K-Kbar
 *     3: eta-eta
 *     4: eta-eta'
 *  4. f0(500) initial state coupling (real part)
 *  5. f0(500) initial state coupling (imaginary part)
 *  6. f0(980) initial state coupling (real part)
 *  7. f0(980) initial state coupling (imaginary part)
 *  8. f0(1370) initial state coupling (real part)
 *  9. f0(1370) initial state coupling (imaginary part)
 * 10. f0(1500) initial state coupling (real part)
 * 11. f0(1500) initial state coupling (imaginary part)
 * 12. f0(1710) initial state coupling (real part)
 * 13. f0(1710) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKMatrixF0::KopfKMatrixF0(const vector<string> &args): UserAmplitude<KopfKMatrixF0>(args) {
	assert(args.size() == 13);
	m_daughters = pair<string, string>(args[0], args[1]);
	channel = atoi(args[2].c_str());
	bf0500_re = AmpParameter(args[3]);
	bf0500_im = AmpParameter(args[4]);
    bf0980_re = AmpParameter(args[5]);
    bf0980_im = AmpParameter(args[6]);
    bf01370_re = AmpParameter(args[7]);
    bf01370_im = AmpParameter(args[8]);
    bf01500_re = AmpParameter(args[9]);
    bf01500_im = AmpParameter(args[10]);
    bf01710_re = AmpParameter(args[11]);
    bf01710_im = AmpParameter(args[12]);
    registerParameter(bf0500_re);
    registerParameter(bf0500_im);
    registerParameter(bf0980_re);
    registerParameter(bf0980_im);
    registerParameter(bf01370_re);
    registerParameter(bf01370_im);
    registerParameter(bf01500_re);
    registerParameter(bf01500_im);
    registerParameter(bf01710_re);
    registerParameter(bf01710_im);
    gf0500 = SVector5(0.74987, -0.01257, 0.27536, -0.15102, 0.36103);
    gf0980 = SVector5(0.06401, 0.00204, 0.77413, 0.50999, 0.13112);
    gf01370 = SVector5(-0.23417, -0.01032, 0.72283, 0.11934, 0.36792);
    gf01500 = SVector5(0.01270, 0.26700, 0.09214, 0.02742, -0.04025);
    gf01710 = SVector5(-0.14242, 0.22780, 0.15981, 0.16272, -0.17397);
    masses = {0.51461, 0.90630, 1.23089, 1.46104, 1.69611};
    couplings = {gf0500, gf0980, gf01370, gf01500, gf01710};
    m1s = {0.1349768, 2*0.1349768, 0.493677, 0.547862, 0.547862};
    m2s = {0.1349768, 2*0.1349768, 0.497611, 0.547862, 0.95778};
    a_bkg = {
        0.03728,
        0.00000, 0.00000,
        -0.01398, 0.00000, 0.02349,
        -0.02203, 0.00000, 0.03101, -0.13769,
        0.01397, 0.00000, -0.04003, -0.06722, -0.28401};
    mat_bkg = SMatrix5Sym(a_bkg.begin(), a_bkg.end());
}


void KopfKMatrixF0::calcUserVars(GDouble** pKin, GDouble* userVars) const {
    TLorentzVector pTemp, pTot;
    /* This allows us to input something like
     * "12 34" for the daughter particle parameters
     * to indicate that the first daughter is the
     * sum of the final state particles indexed
     * 1 and 2 and the second daughter is created
     * from the sum of 3 and 4
     */
    for(unsigned int i=0; i < m_daughters.first.size(); ++i) {
        string num; num += m_daughters.first[i];
        int index = atoi(num.c_str());
        pTemp.SetPxPyPzE(pKin[index][1],
                         pKin[index][2],
                         pKin[index][3],
                         pKin[index][0]);
        pTot += pTemp;
    }
    for(unsigned int i=0; i < m_daughters.second.size(); ++i) {
        string num; num += m_daughters.second[i];
        int index = atoi(num.c_str());
        pTemp.SetPxPyPzE(pKin[index][1],
                         pKin[index][2],
                         pKin[index][3],
                         pKin[index][0]);
        pTot += pTemp;
    }
    // Grab our hypothesis mass
    GDouble m = pTot.M();
    userVars[kM] = m;
    GDouble s = pTot.M2();
    userVars[kS] = s;
    // Calculate K-Matrix
    SMatrix5 mat_K; // Initialized as a 5x5 0-matrix
    SMatrix5 mat_C;
    // Loop over resonances
    for (int i = 0; i < 5; i++) {
        SMatrix5 temp_K;
        SMatrix5 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        temp_K += (mat_bkg * (masses[i] * masses[i] - s));
        temp_K *= KopfKMatrixF0::poleProductRemainder(s, i); 
        // Loop over channels:
        SVector5 B_factor;
        for (int j = 0; j < 5; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        // Loop over channels: 
        for (int k = 0; k < 5; k++) {
            // Loop over channels: 
            for (int l = 0; l < 5; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 5; j++) {
        mat_C(j, j) = KopfKMatrixF0::chew(s, m1s[j], m2s[j]);
    }
    SMatrix5 temp;
    complex<GDouble> product = KopfKMatrixF0::poleProduct(s);
    // Loop over channels
    for (int i = 0; i < 5; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    mat_K *= ((s - 0.0091125) / 1); // Adler zero term
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix5 temp_inv = KopfKMatrixF0::inverse5(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 5; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 5] = temp_inv(channel, i).imag(); // +5 to skip the real parts 
    }
}

complex<GDouble> KopfKMatrixF0::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector5 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bf0500_re, bf0500_im),
        complex<GDouble>(bf0980_re, bf0980_im),
        complex<GDouble>(bf01370_re, bf01370_im),
        complex<GDouble>(bf01500_re, bf01500_im),
        complex<GDouble>(bf01710_re, bf01710_im)
    };
    // Loop over resonances
    for (int i = 0; i < 5; i++) {
        SVector5 temp_P;
        SMatrix5 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        temp_P *= KopfKMatrixF0::poleProductRemainder(s, i);
        SVector5 B_factor;
        // Loop over channels:
        for (int j = 0; j < 5; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        temp_P = temp_P * B_factor;
        vec_P += temp_P;
    }
    SVector5 temp_inv;
    // Loop over channels:
    for (int i = 0; i < 5; i++) {
        temp_inv(i) = complex<GDouble>(userVars[i + 2], userVars[i + 2 + 5]);
    }
    return Dot(temp_inv, vec_P);
}

SMatrix5 KopfKMatrixF0::inverse5(SMatrix5 A) const {
    // A matrix inverse that works for complex values
    SMatrix5 M;
    SMatrix5 I = SMatrixIdentity();
    SMatrix5 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 5; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KopfKMatrixF0::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKMatrixF0::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKMatrixF0::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKMatrixF0::rho(s, m1, m2) / Pi()) * log((KopfKMatrixF0::xi(s, m1, m2) + KopfKMatrixF0::rho(s, m1, m2)) / (KopfKMatrixF0::xi(s, m1, m2) - KopfKMatrixF0::rho(s, m1, m2)));
    tot -= (KopfKMatrixF0::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKMatrixF0::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKMatrixF0::poleProductRemainder(double s, size_t index) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        if (i == index) continue;
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

void KopfKMatrixF0::updatePar(const AmpParameter &par) {}
