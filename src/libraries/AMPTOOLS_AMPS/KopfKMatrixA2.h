#if !defined(KOPFKMATRIXA2)
#define KOPFKMATRIXA2

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

class Kinematics;

class KopfKMatrixA2: public UserAmplitude<KopfKMatrixA2> {
    public:
        KopfKMatrixA2(): UserAmplitude<KopfKMatrixA2>() {}
        KopfKMatrixA2(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re,
            k0im, k1im, k2im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KopfKMatrixA2(){}
    
        string name() const {return "KopfKMatrixA2";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter ba21320_re;
        AmpParameter ba21320_im;
        AmpParameter ba21700_re;
        AmpParameter ba21700_im;

        SVector3 ga21320;
        SVector3 ga21700;
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
