//Jan 18th 2020, Based on model by Adam Szczepaniak & Vincent Mathieu
#if !defined(OMEGAPI_AMPLITUDE)
#define OMEGAPI_AMPLITUDE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"

#ifdef GPU_ACCELERATION
void GPUomegapi_amplitude_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			     int sign, int lambda_gamma, int spin, int parity, int spin_proj, int l,
	 GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction);

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class omegapi_amplitude : public UserAmplitude< omegapi_amplitude >
{

public:
  
  omegapi_amplitude() : UserAmplitude< omegapi_amplitude >() { }
  omegapi_amplitude( const vector< string >& args );
  ~omegapi_amplitude(){}
  
  string name() const { return "omegapi_amplitude"; }
  
complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
//complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  // **********************
  // The following lines are optional and can be used to precalcualte
  // user-defined data that the amplitudes depend on.
  
  // Use this for indexing a user-defined data array and notifying
  // the framework of the number of user-defined variables.
  
  enum UserVars { uv_cosTheta = 0, uv_Phi = 1, uv_cosThetaH = 2, uv_PhiH = 3, uv_prod_angle = 4, uv_Pgamma = 5, uv_dalitz_z = 6, uv_dalitz_sin3theta = 7, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  
  // This function needs to be defined -- see comments and discussion
  // in the .cc file.
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
  
  // This is an optional addition if the calcAmplitude routine
  // can run with only the user-defined data and not the original
  // four-vectors.  It is used to optimize memory usage in GPU
  // based fits.
  bool needsUserVarsOnly() const { return true; }
  // **  end of optional lines **
  
 void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:

  int sign;
  int lambda_gamma;
  int spin;
  int parity;
  int spin_proj;
  int l;
  int nat_sign;
	AmpParameter dalitz_alpha;
	AmpParameter dalitz_beta;
	AmpParameter dalitz_gamma;
	AmpParameter dalitz_delta;

	double polAngle, polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;

};

#endif
