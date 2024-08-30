#if !defined(KOPFKMATRIXA0)
#define KOPFKMATRIXA0

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

class Kinematics;

class KopfKMatrixA0: public UserAmplitude<KopfKMatrixA0> {
    public:
        KopfKMatrixA0(): UserAmplitude<KopfKMatrixA0>() {}
        KopfKMatrixA0(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re,
            k0im, k1im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KopfKMatrixA0(){}
    
        string name() const {return "KopfKMatrixA0";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter ba0980_re;
        AmpParameter ba0980_im;
        AmpParameter ba01450_re;
        AmpParameter ba01450_im;

        SVector2 ga0980;
        SVector2 ga01450;
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
