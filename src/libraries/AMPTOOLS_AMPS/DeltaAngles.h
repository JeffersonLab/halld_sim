#if !defined(DELTAANGLES)
#define DELTAANGLES

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
void GPUDeltaAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			GDouble rho011, GDouble rho031, GDouble rho03m1, 
			GDouble rho111, GDouble rho133, GDouble rho131, GDouble rho13m1,
			GDouble rho231, GDouble rho23m1, GDouble polAngle );
#endif // GPU_ACCELERATION


class Kinematics;

class DeltaAngles : public UserAmplitude< DeltaAngles >
{
    
public:
	
	DeltaAngles() : UserAmplitude< DeltaAngles >() { };
	DeltaAngles( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosTheta, kSinSqTheta, kSin2Theta, kBigPhi, kPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }
	
	string name() const { return "DeltaAngles"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

	bool needsUserVarsOnly() const { return true; }

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
  
	AmpParameter rho011;
	AmpParameter rho031;
	AmpParameter rho03m1;
	
	AmpParameter rho111;
	AmpParameter rho133;
	AmpParameter rho131;
	AmpParameter rho13m1;
	
	AmpParameter rho231;
	AmpParameter rho23m1;

	GDouble polFraction=0.;
	GDouble polAngle=-1;
	TH1D *polFrac_vs_E;

	string lowerVertex;
	string upperVertex;

};

#endif
