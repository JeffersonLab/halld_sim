#if !defined(KOPFKMATRIXF2)
#define KOPFKMATRIXF2

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
typedef SMatrix<complex<GDouble>, 4> SMatrix4;
typedef SMatrix<complex<GDouble>, 4, 4, ROOT::Math::MatRepSym<complex<GDouble>, 4>> SMatrix4Sym;
typedef SVector<complex<GDouble>, 4> SVector4;

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

class Kinematics;

class KopfKMatrixF2: public UserAmplitude<KopfKMatrixF2> {
    public:
        KopfKMatrixF2(): UserAmplitude<KopfKMatrixF2>() {}
        KopfKMatrixF2(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re, k3re,
            k0im, k1im, k2im, k3im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KopfKMatrixF2(){}
    
        string name() const {return "KopfKMatrixF2";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bf21270_re;
        AmpParameter bf21270_im;
        AmpParameter bf21525_re;
        AmpParameter bf21525_im;
        AmpParameter bf21810_re;
        AmpParameter bf21810_im;
        AmpParameter bf21950_re;
        AmpParameter bf21950_im;

        SVector4 gf21270;
        SVector4 gf21525;
        SVector4 gf21810;
        SVector4 gf21950;
        vector<SVector4> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix4 mat_bkg;

        SMatrix4 inverse4(SMatrix4 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
