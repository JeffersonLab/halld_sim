#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <math.h>

#include "TLorentzVector.h"

#include "IUAmpTools/Kinematics.h"

#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/KopfKMatrixF2.h"

/* 
 *  Usage: KopfKMatrixF2 <daughter 1> <daughter 2> <channel> <Re[f₂(1270)]> <Im[f₂(1270)]> ...
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
 *     ╭────┬──────┬────┬────╮
 *     │ 0  │ 1    │ 2  │ 3  │
 *     ├────┼──────┼────┼────┤
 *     │ ππ │ 2π2π │ KK̅ │ ηη │
 *     ╰────┴──────┴────┴────╯
 *  4. f₂(1270) initial state coupling (real part)
 *  5. f₂(1270) initial state coupling (imaginary part)
 *  6. f₂'(1525) initial state coupling (real part)
 *  7. f₂'(1525) initial state coupling (imaginary part)
 *  8. f₂(1810) initial state coupling (real part)
 *  9. f₂(1810) initial state coupling (imaginary part)
 * 10. f₂(1950) initial state coupling (real part)
 * 11. f₂(1950) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

KopfKMatrixF2::KopfKMatrixF2(const vector<string> &args): UserAmplitude<KopfKMatrixF2>(args) {
	assert(args.size() == 11);
	m_daughters = pair<string, string>(args[0], args[1]);
	channel = atoi(args[2].c_str());
	bf21270_re = AmpParameter(args[3]);
	bf21270_im = AmpParameter(args[4]);
    bf21525_re = AmpParameter(args[5]);
    bf21525_im = AmpParameter(args[6]);
    bf21810_re = AmpParameter(args[7]);
    bf21810_im = AmpParameter(args[8]);
    bf21950_re = AmpParameter(args[9]);
    bf21950_im = AmpParameter(args[10]);
    registerParameter(bf21270_re);
    registerParameter(bf21270_im);
    registerParameter(bf21525_re);
    registerParameter(bf21525_im);
    registerParameter(bf21810_re);
    registerParameter(bf21810_im);
    registerParameter(bf21950_re);
    registerParameter(bf21950_im);
    gf21270 = SVector4(0.40033, 0.15479, -0.08900, -0.00113);
    gf21525 = SVector4(0.01820, 0.17300, 0.32393, 0.15256);
    gf21810 = SVector4(-0.06709, 0.22941, -0.43133, 0.23721);
    gf21950 = SVector4(-0.49924, 0.19295, 0.27975, -0.03987);
    masses = {1.15299, 1.48359, 1.72923, 1.96700};
    couplings = {gf21270, gf21525, gf21810, gf21950};
    m1s = {0.1349768, 2*0.1349768, 0.493677, 0.547862};
    m2s = {0.1349768, 2*0.1349768, 0.497611, 0.547862};
    a_bkg = {
        -0.04319,
        0.00000, 0.00000,
        0.00984, 0.00000, -0.07344,
        0.01028, 0.00000, 0.05533, -0.05183};
    mat_bkg = SMatrix4Sym(a_bkg.begin(), a_bkg.end());
}

void KopfKMatrixF2::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix4 mat_K; // Initialized as a 4x4 0-matrix
    SMatrix4 mat_C;
    // Loop over resonances
    for (int i = 0; i < 4; i++) {
        SMatrix4 temp_K;
        SMatrix4 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        temp_K += (mat_bkg * (masses[i] * masses[i] - s));
        temp_K *= KopfKMatrixF2::poleProductRemainder(s, i);
        // Loop over channels: 
        SVector4 B_factor;
        for (int j = 0; j < 4; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        // Loop over channels: 
        for (int k = 0; k < 4; k++) {
            // Loop over channels: 
            for (int l = 0; l < 4; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 4; j++) {
        mat_C(j, j) = KopfKMatrixF2::chew(s, m1s[j], m2s[j]);
    }
    SMatrix4 temp;
    complex<GDouble> product = KopfKMatrixF2::poleProduct(s);
    // Loop over channels
    for (int i = 0; i < 4; ++i) {
        temp(i, i) = product;
    }
    mat_K *= mat_C;
    temp += mat_K;
    // Now temp is (I + KC)
    SMatrix4 temp_inv = KopfKMatrixF2::inverse4(temp); // (I + KC)^{-1}
    // Now we cache the results
    for (int i = 0; i < 4; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 4] = temp_inv(channel, i).imag(); // +4 to skip the real parts 
    }
}



complex<GDouble> KopfKMatrixF2::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector4 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bf21270_re, bf21270_im),
        complex<GDouble>(bf21525_re, bf21525_im),
        complex<GDouble>(bf21810_re, bf21810_im),
        complex<GDouble>(bf21950_re, bf21950_im)
    };
    // Loop over resonances
    for (int i = 0; i < 4; i++) {
        SVector4 temp_P;
        SMatrix4 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        temp_P *= KopfKMatrixF2::poleProductRemainder(s, i);
        // Loop over channels:
        SVector4 B_factor;
        for (int j = 0; j < 4; j++) {
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
    SVector4 temp_inv;
    // Loop over channels:
    for (int i = 0; i < 4; i++) {
        temp_inv(i) = complex<GDouble>(userVars[i + 2], userVars[i + 2 + 4]);
    }
    return Dot(temp_inv, vec_P);
}

SMatrix4 KopfKMatrixF2::inverse4(SMatrix4 A) const {
    // A matrix inverse that works for complex values
    SMatrix4 M;
    SMatrix4 I = SMatrixIdentity();
    SMatrix4 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 4; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KopfKMatrixF2::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KopfKMatrixF2::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KopfKMatrixF2::chew(double s, double m1, double m2) const {
    // The Chew-Mandelstam matrix as described by Appendix B in <https://doi.org/10.1103/PhysRevD.91.054008>
    complex<GDouble> tot = 0;
    tot += (KopfKMatrixF2::rho(s, m1, m2) / Pi()) * log((KopfKMatrixF2::xi(s, m1, m2) + KopfKMatrixF2::rho(s, m1, m2)) / (KopfKMatrixF2::xi(s, m1, m2) - KopfKMatrixF2::rho(s, m1, m2)));
    tot -= (KopfKMatrixF2::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

complex<GDouble> KopfKMatrixF2::poleProduct(double s) const {
    // We multiply the numerator and denominator by the product of poles to remove
    // division by zero errors when the measured mass is computationally close to
    // the pole mass.
    complex<GDouble> tot = 1;
    for (size_t i = 0; i < masses.size(); ++i) {
        tot *= complex<GDouble>(masses[i] * masses[i] - s);
    }
    return tot;
}

complex<GDouble> KopfKMatrixF2::poleProductRemainder(double s, size_t index) const {
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

void KopfKMatrixF2::updatePar(const AmpParameter &par) {}
