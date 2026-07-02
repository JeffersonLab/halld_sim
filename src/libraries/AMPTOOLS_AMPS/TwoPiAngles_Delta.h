#if !defined(TWOPIANGLES_DELTA)
#define TWOPIANGLES_DELTA

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
GPUTwoPiAngles_Delta_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
			 GDouble rho000, GDouble rho100, GDouble rho1m10,
			 GDouble rho111, GDouble delta_rho133, GDouble rho101,
			 GDouble rho1m11, GDouble rho102, GDouble rho1m12,
			 GDouble polAngle,GDouble delta_rho111 );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class TwoPiAngles_Delta : public UserAmplitude< TwoPiAngles_Delta >
{
    
public:
	
	TwoPiAngles_Delta() : UserAmplitude< TwoPiAngles_Delta >() { };
	TwoPiAngles_Delta( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosTheta, kSinSqTheta, kSin2Theta,
			kBigPhi, kPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }

	string name() const { return "TwoPiAngles_Delta"; }
    
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
	
  AmpParameter rho111;
  AmpParameter delta_rho133;
  AmpParameter rho101;
  AmpParameter rho1m11;

  AmpParameter rho102;
  AmpParameter rho1m12;
  AmpParameter rho001;
  
  AmpParameter polAngle;
  AmpParameter delta_rho111;
  bool constraint;

  double polFraction;
  //TH1D *polFrac_vs_E;
  string frame;

};

#endif
