#if !defined(TWOPIANGLES_DELTA_FACTORIZED)
#define TWOPIANGLES_DELTA_FACTORIZED

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
GPUTwoPiAngles_Delta_factorized_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
  GDouble rho000,
  GDouble rho100,
  GDouble rho1m10,
  GDouble rho001,
  GDouble rho111,
  GDouble rho101,
  GDouble rho1m11,
  GDouble rho102,
  GDouble rho1m12,
  GDouble delta_rho011,
  GDouble delta_rho031,
  GDouble delta_rho03m1
);

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class TwoPiAngles_Delta_factorized : public UserAmplitude< TwoPiAngles_Delta_factorized >
{
    
public:
	
  TwoPiAngles_Delta_factorized() : UserAmplitude< TwoPiAngles_Delta_factorized >() { };
	TwoPiAngles_Delta_factorized( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosTheta_pim, kSinSqTheta_pim, kSin2Theta_pim, kPhi_pim, kCosTheta_proton, kSinSqTheta_proton, kSin2Theta_proton, kPhi_proton,
			kBigPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }

	string name() const { return "TwoPiAngles_Delta_factorized"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

	// we can calcualte everythign we need from userVars block so allow
	// the framework to purge the four-vectors
	bool needsUserVarsOnly() const { return true; }
	
#ifdef GPU_ACCELERATION
  
	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:

AmpParameter rho000;
AmpParameter rho100;
AmpParameter rho1m10;

AmpParameter rho001;
AmpParameter rho111;
AmpParameter rho101;
AmpParameter rho1m11;

AmpParameter rho102;
AmpParameter rho1m12;
AmpParameter delta_rho011;
AmpParameter delta_rho031;
AmpParameter delta_rho03m1;

AmpParameter polAngle;
  
  double polFraction;
  //TH1D *polFrac_vs_E;
  string frame;
};

#endif
