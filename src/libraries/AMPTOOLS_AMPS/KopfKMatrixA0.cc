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
#include "AMPTOOLS_AMPS/KopfKMatrixA0.h"

/* 
 *  Usage: KopfKMatrixA0 <daughter 1> <daughter 2> <channel> <Re[a₀(980)]> <Im[a₀(980)]> ...
 * 
 *  Parameters:
 *  1. <daughter 1>: a single digit number or pair of digits corresponding to the location of
 *     the first decay product in the kinematic arrays. Note that this array typically starts
 *     with the beam in location 0, so these values could be considered to be 1-indexed. If
 *     the number is more than one digit long, each single digit will represent a particle
 *     which, when added together, form the first decay product (e.g. 23 would sum particles
 *     2 and 3 to make a new state).
 *  2. <daughter 2>: Same as (1) but for the second decay product.
 *  3. <channel>: A single digit corresponding to the final state channel (see table).
 *     ╭────┬────╮
 *     │ 0  │ 1  │
 *     ├────┼────┤
 *     │ ηπ │ KK̅ │
 *     ╰────┴────╯
 *  4. a₀(980) initial state coupling (real part)
 *  5. a₀(980) initial state coupling (imaginary part)
 *  6. a₀(1450) initial state coupling (real part)
 *  7. a₀(1450) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKmatrixA0::KopfKmatrixA0(const vector<string> &args): UserAmplitude<KopfKmatrixA0>(args) {
	assert(args.size() == 7);
	m_daughters = pair<string, string>(args[0], args[1]);
	channel = atoi(args[2].c_str());
	ba0980_re = AmpParameter(args[3]);
	ba0980_im = AmpParameter(args[4]);
    ba01450_re = AmpParameter(args[5]);
    ba01450_im = AmpParameter(args[6]);
    registerParameter(ba0980_re);
    registerParameter(ba0980_im);
    registerParameter(ba01450_re);
    registerParameter(ba01450_im);
    ga0980 = SVector2(0.43215, -0.28825);
    ga01450 = SVector2(0.19000, 0.43372);
    masses = {0.95395, 1.26767};
    couplings = {ga0980, ga01450};
    m1s = {0.1349768, 0.493677};
    m2s = {0.547862, 0.497611};
    a_bkg = {
        0.00000,
        0.00000, 0.00000};
    mat_bkg = SMatrix2Sym(a_bkg.begin(), a_bkg.end());
}

void KopfKmatrixA0::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix2 mat_K; // Initialized as a 4x4 0-matrix
    SMatrix2 mat_C;
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SMatrix2 temp_K;
        SMatrix2 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        temp_K += (mat_bkg * (masses[i] * masses[i] - s));
        temp_K *= KopfKmatrixA0::poleProductRemainder(s, i);
        // Loop over channels: 
        SVector2 B_factor;
        for (int j = 0; j < 2; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        // Loop over channels: 
        for (int k = 0; k < 2; k++) {
            // Loop over channels: 
            for (int l = 0; l < 2; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 2; j++) {
        mat_C(j, j) = KopfKmatrixA0::chew(s, m1s[j], m2s[j]);
    }
    SMatrix2 temp;
    complex<GDouble> product = KopfKmatrixA0::poleProduct(s);
    for (int i = 0; i < 2; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix2 temp_inv = KopfKmatrixA0::inverse2(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 2; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 2] = temp_inv(channel, i).imag(); // +2 to skip the real parts 
    }
}



complex<GDouble> KopfKmatrixA0::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector2 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(ba0980_re, ba0980_im),
        complex<GDouble>(ba01450_re, ba01450_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector2 temp_P;
        SMatrix2 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        temp_P *= KopfKmatrixA0::poleProductRemainder(s, i);
        SVector2 B_factor;
        // Loop over channels:
        for (int j = 0; j < 2; j++) {
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
    SVector2 temp_inv;
    // Loop over channels:
    for (int i = 0; i < 2; i++) {
        temp_inv(i) = complex<GDouble>(userVars[i + 2], userVars[i + 2 + 2]);
    }
    return Dot(temp_inv, vec_P);
}

SMatrix2 KopfKmatrixA0::inverse2(SMatrix2 A) const {
    // A matrix inverse that works for complex values
    SMatrix2 M;
    SMatrix2 I = SMatrixIdentity();
    SMatrix2 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 2; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KopfKmatrixA0::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKmatrixA0::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKmatrixA0::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKmatrixA0::rho(s, m1, m2) / Pi()) * log((KopfKmatrixA0::xi(s, m1, m2) + KopfKmatrixA0::rho(s, m1, m2)) / (KopfKmatrixA0::xi(s, m1, m2) - KopfKmatrixA0::rho(s, m1, m2)));
    tot -= (KopfKmatrixA0::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKmatrixA0::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKmatrixA0::poleProductRemainder(double s, size_t index) const {
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

void KopfKmatrixA0::updatePar(const AmpParameter &par) {}
