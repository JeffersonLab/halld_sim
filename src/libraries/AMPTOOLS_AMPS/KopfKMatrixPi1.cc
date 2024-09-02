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
#include "AMPTOOLS_AMPS/KopfKMatrixPi1.h"

/* 
 *  Usage: KopfKMatrixPi1 <daughter 1> <daughter 2> <channel> <Re[π1(1600)]> <Im[π1(1600)]> ...
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
 *     1: pi-eta'
 *  4. π1(1600) initial state coupling (real part)
 *  5. π1(1600) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKMatrixPi1::KopfKMatrixPi1(const vector<string> &args): UserAmplitude<KopfKMatrixPi1>(args) {
	assert(args.size() == 5);
    m_daughters = pair<string, string>(args[0], args[1]);
    channel = atoi(args[2].c_str());
    bpi11600_re = AmpParameter(args[3]);
    bpi11600_im = AmpParameter(args[4]);
    registerParameter(bpi11600_re);
    registerParameter(bpi11600_im);
    gpi11600 = SVector2(0.80564, 1.04695);
    masses = {1.38552};
    couplings = {gpi11600};
    m1s = {0.1349768, 0.1349768};
    m2s = {0.547862, 0.95778};
    a_bkg = {
        1.05000,
        0.15163, -0.24611};
    mat_bkg = SMatrix2Sym(a_bkg.begin(), a_bkg.end());
}

void KopfKMatrixPi1::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
        temp_K *= KopfKMatrixPi1::poleProductRemainder(s, i);
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
        mat_C(j, j) = KopfKMatrixPi1::chew(s, m1s[j], m2s[j]);
    }
    SMatrix2 temp;
    complex<GDouble> product = KopfKMatrixPi1::poleProduct(s);
    // Loop over channels
    for (int i = 0; i < 2; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix2 temp_inv = KopfKMatrixPi1::inverse2(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 2; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 2] = temp_inv(channel, i).imag(); // +2 to skip the real parts 
    }
}



complex<GDouble> KopfKMatrixPi1::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector2 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bpi11600_re, bpi11600_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector2 temp_P;
        SMatrix2 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        temp_P *= KopfKMatrixPi1::poleProductRemainder(s, i);
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

SMatrix2 KopfKMatrixPi1::inverse2(SMatrix2 A) const {
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

complex<GDouble> KopfKMatrixPi1::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKMatrixPi1::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKMatrixPi1::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKMatrixPi1::rho(s, m1, m2) / Pi()) * log((KopfKMatrixPi1::xi(s, m1, m2) + KopfKMatrixPi1::rho(s, m1, m2)) / (KopfKMatrixPi1::xi(s, m1, m2) - KopfKMatrixPi1::rho(s, m1, m2)));
    tot -= (KopfKMatrixPi1::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKMatrixPi1::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKMatrixPi1::poleProductRemainder(double s, size_t index) const {
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

void KopfKMatrixPi1::updatePar(const AmpParameter &par) {}
