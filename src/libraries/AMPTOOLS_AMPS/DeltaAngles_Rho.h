#if !defined(DELTAANGLES_RHO)
#define DELTAANGLES_RHO

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
void GPUDeltaAngles_Rho_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			GDouble delta_rho011, GDouble delta_rho031, GDouble delta_rho03m1, 
			GDouble delta_rho111, GDouble delta_rho133, GDouble delta_rho131, GDouble delta_rho13m1,
			GDouble delta_rho231, GDouble delta_rho23m1, GDouble polAngle );
#endif // GPU_ACCELERATION

class Kinematics;

class DeltaAngles_Rho : public UserAmplitude< DeltaAngles_Rho >
{
    
public:
	
	DeltaAngles_Rho() : UserAmplitude< DeltaAngles_Rho >() { };
	DeltaAngles_Rho( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosTheta, kSinSqTheta, kSin2Theta, kBigPhi, kPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }

	string name() const { return "DeltaAngles_Rho"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars  ) const;
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

	bool needsUserVarsOnly() const { return true; }

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
  
	AmpParameter delta_rho011;
	AmpParameter delta_rho031;
	AmpParameter delta_rho03m1;
	
	AmpParameter delta_rho111;
	AmpParameter delta_rho133;
	AmpParameter delta_rho131;
	AmpParameter delta_rho13m1;
	
	AmpParameter delta_rho231;
	AmpParameter delta_rho23m1;

	GDouble polFraction=0.;
	GDouble polAngle=-1;
	TH1D *polFrac_vs_E;
	string frame;

};

#endif