#if !defined(KOPFKMATRIXRHO)
#define KOPFKMATRIXRHO

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
typedef SMatrix<complex<GDouble>, 3> SMatrix3;
typedef SMatrix<complex<GDouble>, 3, 3, ROOT::Math::MatRepSym<complex<GDouble>, 3>> SMatrix3Sym;
typedef SVector<complex<GDouble>, 3> SVector3;

/* 
 *  Usage: KopfKMatrixRho <daughter 1> <daughter 2> <channel> <Re[rho(770)]> <Im[rho(770)]> ...
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
 *  4. rho(770) initial state coupling (real part)
 *  5. rho(770) initial state coupling (imaginary part)
 *  6. rho(1700) initial state coupling (real part)
 *  7. rho(1700) initial state coupling (imaginary part)
 *
 *  See <https://doi.org/10.1140/epjc/s10052-021-09821-2> for more details.
 */

class Kinematics;

class KopfKMatrixRho: public UserAmplitude<KopfKMatrixRho> {
    public:
        KopfKMatrixRho(): UserAmplitude<KopfKMatrixRho>() {}
        KopfKMatrixRho(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re,
            k0im, k1im, k2im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KopfKMatrixRho(){}
    
        string name() const {return "KopfKMatrixRho";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter brho770_re;
        AmpParameter brho770_im;
        AmpParameter brho1700_re;
        AmpParameter brho1700_im;

        SVector3 grho770;
        SVector3 grho1700;
        vector<SVector3> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix3 mat_bkg;

        SMatrix3 inverse3(SMatrix3 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
