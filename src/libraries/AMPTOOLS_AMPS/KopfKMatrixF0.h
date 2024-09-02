#if !defined(KOPFKMATRIXF0)
#define KOPFKMATRIXF0

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include "Math/SVector.h"
#include "Math/SMatrix.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;
using namespace ROOT::Math;
typedef SMatrix<complex<GDouble>, 5> SMatrix5;
typedef SMatrix<complex<GDouble>, 5, 5, ROOT::Math::MatRepSym<complex<GDouble>, 5>> SMatrix5Sym;
typedef SVector<complex<GDouble>, 5> SVector5;

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

class Kinematics;

class KopfKMatrixF0: public UserAmplitude<KopfKMatrixF0> {
    public:
        KopfKMatrixF0(): UserAmplitude<KopfKMatrixF0>() {}
        KopfKMatrixF0(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re, k3re, k4re,
            k0im, k1im, k2im, k3im, k4im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}
        
        ~KopfKMatrixF0(){}
    
        string name() const {return "KopfKMatrixF0";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bf0500_re;
        AmpParameter bf0500_im;
        AmpParameter bf0980_re;
        AmpParameter bf0980_im;
        AmpParameter bf01370_re;
        AmpParameter bf01370_im;
        AmpParameter bf01500_re;
        AmpParameter bf01500_im;
        AmpParameter bf01710_re;
        AmpParameter bf01710_im;

        SVector5 gf0500;
        SVector5 gf0980;
        SVector5 gf01370;
        SVector5 gf01500;
        SVector5 gf01710;
        vector<SVector5> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix5 mat_bkg;

        SMatrix5 inverse5(SMatrix5 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
