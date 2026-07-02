#if !defined(TWOPIANGLES_DELTA_DOUBLESDMES_UNPOL)
#define TWOPIANGLES_DELTA_DOUBLESDMES_UNPOL

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
GPUTwoPiAngles_Delta_DoubleSDMEs_unpol_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
GDouble r00_33_0, //GDouble r00_11_0
GDouble r11_33_0, GDouble r11_11_0,
GDouble r1m1_33_0, GDouble r1m1_11_0,
GDouble r10_33_0, GDouble r10_11_0,
GDouble r00_31_0, GDouble r00_3m1_0,
GDouble r11_31_0, GDouble r11_3m1_0,
GDouble r1m1_31_0, GDouble r1m1_3m1_0,
GDouble r10_31_0, GDouble r10_3m1_0,
GDouble rt1m1_31_0, GDouble rt1m1_3m1_0,
GDouble rt10_31_0,  GDouble rt10_3m1_0

);

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class TwoPiAngles_Delta_DoubleSDMEs_unpol : public UserAmplitude< TwoPiAngles_Delta_DoubleSDMEs_unpol >
{
    
public:
	
	TwoPiAngles_Delta_DoubleSDMEs_unpol() : UserAmplitude< TwoPiAngles_Delta_DoubleSDMEs_unpol >() { };
	TwoPiAngles_Delta_DoubleSDMEs_unpol( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosTheta_pim, kSinSqTheta_pim, kSin2Theta_pim, kPhi_pim, kCosTheta_proton, kSinSqTheta_proton, kSin2Theta_proton, kPhi_proton,
			kBigPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }

	string name() const { return "TwoPiAngles_Delta_DoubleSDMEs_unpol"; }
    
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

  AmpParameter r00_33_0;
  // AmpParameter r00_11_0;

  AmpParameter r11_33_0;
  AmpParameter r11_11_0;

  AmpParameter r1m1_33_0;
  AmpParameter r1m1_11_0;

  AmpParameter r10_33_0;
  AmpParameter r10_11_0;

  // \bar W^0 (31 / 3-1)

  AmpParameter r00_31_0;
  AmpParameter r00_3m1_0;

  AmpParameter r11_31_0;
  AmpParameter r11_3m1_0;

  AmpParameter r1m1_31_0;
  AmpParameter r1m1_3m1_0;

  AmpParameter r10_31_0;
  AmpParameter r10_3m1_0;

  // \tilde W^0

  AmpParameter rt1m1_31_0;
  AmpParameter rt1m1_3m1_0;

  AmpParameter rt10_31_0;
  AmpParameter rt10_3m1_0;

  AmpParameter polAngle;
  
  double polFraction;
  //TH1D *polFrac_vs_E;
  string frame;
  string AMO;
};

#endif
