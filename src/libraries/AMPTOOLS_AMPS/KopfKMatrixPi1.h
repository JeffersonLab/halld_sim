#if !defined(KOPFKMATRIXPI1)
#define KOPFKMATRIXPI1

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
typedef SMatrix<complex<GDouble>, 2> SMatrix2;
typedef SMatrix<complex<GDouble>, 2, 2, ROOT::Math::MatRepSym<complex<GDouble>, 2>> SMatrix2Sym;
typedef SVector<complex<GDouble>, 2> SVector2;

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

class Kinematics;

class KopfKMatrixPi1: public UserAmplitude<KopfKMatrixPi1> {
    public:
        KopfKMatrixPi1(): UserAmplitude<KopfKMatrixPi1>() {}
        KopfKMatrixPi1(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re,
            k0im, k1im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KopfKMatrixPi1(){}
    
        string name() const {return "KopfKMatrixPi1";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bpi11600_re;
        AmpParameter bpi11600_im;

        SVector2 gpi11600;
        vector<SVector2> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix2 mat_bkg;

        SMatrix2 inverse2(SMatrix2 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
