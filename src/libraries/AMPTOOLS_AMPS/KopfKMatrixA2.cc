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
#include "AMPTOOLS_AMPS/KopfKMatrixA2.h"

/* 
 *  Usage: KopfKMatrixA2 <daughter 1> <daughter 2> <channel> <Re[a2(1320)]> <Im[a2(1320)]> ...
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
 *     0: pi-eta
 *     1: K-Kbar
 *     2: pi-eta'
 *  4. a2(1320) initial state coupling (real part)
 *  5. a2(1320) initial state coupling (imaginary part)
 *  6. a2(1700) initial state coupling (real part)
 *  7. a2(1700) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKMatrixA2::KopfKMatrixA2(const vector<string> &args): UserAmplitude<KopfKMatrixA2>(args) {
	assert(args.size() == 7);
	m_daughters = pair<string, string>(args[0], args[1]);
	channel = atoi(args[2].c_str());
	ba21320_re = AmpParameter(args[3]);
	ba21320_im = AmpParameter(args[4]);
    ba21700_re = AmpParameter(args[5]);
    ba21700_im = AmpParameter(args[6]);
    registerParameter(ba21320_re);
    registerParameter(ba21320_im);
    registerParameter(ba21700_re);
    registerParameter(ba21700_im);
    ga21320 = SVector3(0.30073, 0.21426, -0.09162);
    ga21700 = SVector3(0.68567, 0.12543, 0.00184);
    masses = {1.30080, 1.75351};
    couplings = {ga21320, ga21700};
    m1s = {0.1349768, 0.493677, 0.1349768};
    m2s = {0.547862, 0.497611, 0.95778};
    a_bkg = {
        -0.40184,
        0.00033, -0.21416,
        -0.08707, -0.06193, -0.17435};
    mat_bkg = SMatrix3Sym(a_bkg.begin(), a_bkg.end());
}

void KopfKMatrixA2::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix3 mat_K; // Initialized as a 3x3 0-matrix
    SMatrix3 mat_C;
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SMatrix3 temp_K;
        SMatrix3 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        temp_K += (mat_bkg * (masses[i] * masses[i] - s));
        temp_K *= KopfKMatrixA2::poleProductRemainder(s, i);
        // Loop over channels: 
        SVector3 B_factor;
        for (int j = 0; j < 3; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        // Loop over channels: 
        for (int k = 0; k < 3; k++) {
            // Loop over channels: 
            for (int l = 0; l < 3; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 3; j++) {
        mat_C(j, j) = KopfKMatrixA2::chew(s, m1s[j], m2s[j]);
    }
    SMatrix3 temp;
    complex<GDouble> product = KopfKMatrixA2::poleProduct(s);
    for (int i = 0; i < 3; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix3 temp_inv = KopfKMatrixA2::inverse3(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 3; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 3] = temp_inv(channel, i).imag(); // +3 to skip the real parts 
    }
}



complex<GDouble> KopfKMatrixA2::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector3 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(ba21320_re, ba21320_im),
        complex<GDouble>(ba21700_re, ba21700_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector3 temp_P;
        SMatrix3 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        temp_P *= KopfKMatrixA2::poleProductRemainder(s, i);
        SVector3 B_factor;
        // Loop over channels:
        for (int j = 0; j < 3; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        temp_P = temp_P * B_factor;
        vec_P += temp_P;
    }
    SVector3 temp_inv;
    // Loop over channels:
    for (int i = 0; i < 3; i++) {
        temp_inv(i) = complex<GDouble>(userVars[i + 2], userVars[i + 2 + 3]);
    }
    return Dot(temp_inv, vec_P);
}

SMatrix3 KopfKMatrixA2::inverse3(SMatrix3 A) const {
    // A matrix inverse that works for complex values
    SMatrix3 M;
    SMatrix3 I = SMatrixIdentity();
    SMatrix3 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 3; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KopfKMatrixA2::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKMatrixA2::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKMatrixA2::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKMatrixA2::rho(s, m1, m2) / Pi()) * log((KopfKMatrixA2::xi(s, m1, m2) + KopfKMatrixA2::rho(s, m1, m2)) / (KopfKMatrixA2::xi(s, m1, m2) - KopfKMatrixA2::rho(s, m1, m2)));
    tot -= (KopfKMatrixA2::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKMatrixA2::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKMatrixA2::poleProductRemainder(double s, size_t index) const {
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

void KopfKMatrixA2::updatePar(const AmpParameter &par) {}
