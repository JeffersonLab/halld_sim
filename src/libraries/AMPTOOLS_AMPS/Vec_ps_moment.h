#if !defined(VEC_PS_MOMENT)
#define VEC_PS_MOMENT

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

#ifdef GPU_ACCELERATION
void
GPUVec_ps_moment_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, moment* moments, int numberOfMoments );
#endif

class Kinematics;

struct moment {
    string name;
    AmpParameter H;
    int alpha;
    int Jv;
    int Lambda;
    int J;
    int M;
};

// An AmpTools class for describing the polarized moments for R-> Vector Pseudoscalar
// with a polarized photon beam, must have m >= 0
class Vec_ps_moment : public UserAmplitude< Vec_ps_moment >
{

public:

    Vec_ps_moment() : UserAmplitude< Vec_ps_moment >() { };
    Vec_ps_moment( const vector< string >& args );

    string name() const { return "Vec_ps_moment"; }

    complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

    enum UserVars { kBeamPolFraction = 0, kCosTheta, kPhi, kCosThetaH, kPhiH, kProdAngle, kNumUserVars };
    unsigned int numUserVars() const { return kNumUserVars; }

    void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

    // we can calcualte everything we need from userVars block so allow
    // the framework to purge the four-vectors
    bool needsUserVarsOnly() const { return true; }

#ifdef GPU_ACCELERATION

    void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

        bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:

    bool m_3pi;    
    int m_numberOfMoments;
    int m_nonMomentArgs;
    vector<moment> m_moments;
    vector<AmpParameter> m_H;    

    double m_polAngle;
    double m_polFraction;
    bool m_polInTree;

    TH1D *m_polFrac_vs_E;
};

#endif
