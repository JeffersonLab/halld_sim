#if !defined(ZLM)
#define ZLM

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUZlm_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
      int j, int m, int r, int s );
#endif // GPU_ACCELERATION


using std::complex;
using namespace std;

// A class for describing the angular portion of the decay R->12
// in the reflectivity basis for photon beams: m can be negative!
// particles 1 and 2 are pseudoscalars
//
// j,m are the total and z projection of the spin of R
// r=+/-1 indicates real/imaginary part of Zlm
// s=+/-1 multiplies with sqrt(1+/- P_gamma)

class Kinematics;

class Zlm : public UserAmplitude< Zlm >
{

   public:

      Zlm() : UserAmplitude< Zlm >() { };
      Zlm( const vector< string >& args );

      enum UserVars { kPgamma = 0, kCosTheta, kPhi, kBigPhi, kNumUserVars };
      unsigned int numUserVars() const { return kNumUserVars; }

      string name() const { return "Zlm"; }

      complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
      void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

      // we can calcualte everything we need from userVars block so allow
      // the framework to purge the four-vectors
      bool needsUserVarsOnly() const { return true; }

      // the user variables above are the same for all instances of this amplitude
      bool areUserVarsStatic() const { return true; }

#ifdef GPU_ACCELERATION

      void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

      bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

   private:

      int m_j;
      int m_m;
      int m_r;
      int m_s;

      double m_polAngle;
      double m_polFraction;
      bool m_polInTree;

      TH1D* m_polFrac_vs_E;
};

#endif
