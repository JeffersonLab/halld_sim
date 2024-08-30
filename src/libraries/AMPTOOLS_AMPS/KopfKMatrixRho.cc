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
#include "AMPTOOLS_AMPS/KopfKMatrixRho.h"

/* 
 *  Usage: KopfKMatrixRho <daughter 1> <daughter 2> <channel> <Re[ρ(770)]> <Im[ρ(770)]> ...
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
 *     ╭────┬──────┬────╮
 *     │ 0  │ 1    │ 2  │
 *     ├────┼──────┼────┤
 *     │ ππ │ 2π2π │ KK̅ │
 *     ╰────┴──────┴────╯
 *  4. ρ(770) initial state coupling (real part)
 *  5. ρ(770) initial state coupling (imaginary part)
 *  6. ρ(1700) initial state coupling (real part)
 *  7. ρ(1700) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKMatrixRho::KopfKMatrixRho(const vector<string> &args): UserAmplitude<KopfKMatrixRho>(args) {
	assert(args.size() == 7);
	m_daughters = pair<string, string>(args[0], args[1]);
	channel = atoi(args[2].c_str());
	brho770_re = AmpParameter(args[3]);
	brho770_im = AmpParameter(args[4]);
    brho1700_re = AmpParameter(args[5]);
    brho1700_im = AmpParameter(args[6]);
    registerParameter(brho770_re);
    registerParameter(brho770_im);
    registerParameter(brho1700_re);
    registerParameter(brho1700_im);
    grho770 = SVector3(0.28023, 0.01806, 0.06501);
    grho1700 = SVector3(0.16318, 0.53879, 0.00495);
    masses = {0.71093, 1.58660};
    couplings = {grho770, grho1700};
    m1s = {0.1349768, 2*0.1349768, 0.493677};
    m2s = {0.1349768, 2*0.1349768, 0.497611};
    a_bkg = {
        -0.06948,
        0.00000, 0.00000,
        0.07958, 0.00000, -0.60000};
    mat_bkg = SMatrix3Sym(a_bkg.begin(), a_bkg.end());
}

void KopfKMatrixRho::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
        temp_K *= KopfKMatrixRho::poleProductRemainder(s, i);
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
        mat_C(j, j) = KopfKMatrixRho::chew(s, m1s[j], m2s[j]);
    }
    SMatrix3 temp;
    complex<GDouble> product = KopfKMatrixRho::poleProduct(s);
    // Loop over channels
    for (int i = 0; i < 3; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix3 temp_inv = KopfKMatrixRho::inverse3(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 3; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 3] = temp_inv(channel, i).imag(); // +3 to skip the real parts 
    }
}



complex<GDouble> KopfKMatrixRho::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector3 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(brho770_re, brho770_im),
        complex<GDouble>(brho1700_re, brho1700_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector3 temp_P;
        SMatrix3 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        temp_P *= KopfKMatrixRho::poleProductRemainder(s, i);
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

SMatrix3 KopfKMatrixRho::inverse3(SMatrix3 A) const {
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

complex<GDouble> KopfKMatrixRho::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKMatrixRho::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKMatrixRho::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKMatrixRho::rho(s, m1, m2) / Pi()) * log((KopfKMatrixRho::xi(s, m1, m2) + KopfKMatrixRho::rho(s, m1, m2)) / (KopfKMatrixRho::xi(s, m1, m2) - KopfKMatrixRho::rho(s, m1, m2)));
    tot -= (KopfKMatrixRho::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKMatrixRho::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKMatrixRho::poleProductRemainder(double s, size_t index) const {
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

void KopfKMatrixRho::updatePar(const AmpParameter &par) {}
