#if !defined(VECPSMOMENT)
#define VECPSMOMENT

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
GPUVecPSMoment_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		     GDouble* H, int* indices, int nMoments );
#endif // GPU_ACCELERATION


using std::complex;
using namespace std;

// A class for describing the polarized moments for R-> Vector Pseudoscalar
// with a polarized photon beam, must have m >= 0
class Kinematics;

class VecPSMoment : public UserAmplitude< VecPSMoment >
{

   public:

      VecPSMoment() : UserAmplitude< VecPSMoment >() { };
      VecPSMoment( const vector< string >& args );

      enum UserVars { kPgamma = 0, kCosTheta, kPhi, kCosThetaH, kPhiH, kBigPhi, kNumUserVars };
      unsigned int numUserVars() const { return kNumUserVars; }

      string name() const { return "VecPSMoment"; }

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

      bool m_3pi;
      int m_maxJ;
      int m_nMoments;
      vector<AmpParameter> m_H;
      vector<int> m_indices;

      double m_polAngle;
      double m_polFraction;
      bool m_polInTree;

      TH1D* m_polFrac_vs_E;
};

#endif
