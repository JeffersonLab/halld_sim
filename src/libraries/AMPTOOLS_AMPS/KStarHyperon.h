#if !defined(KSTARHYPERON)
#define KSTARHYPERON

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
// GPUKStarHyperon_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
// 	    GDouble alpha, GDouble Sigma, GDouble Ox, GDouble P, GDouble T, GDouble Oz,
// 	    GDouble polAngle );
GPUKStarHyperon_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
	GDouble alpha, GDouble rho111, GDouble rho001, GDouble Ox, GDouble P, GDouble T, GDouble Oz,
	GDouble polAngle );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class KStarHyperon : public UserAmplitude< KStarHyperon >
{
    
public:
	
	KStarHyperon() : UserAmplitude< KStarHyperon >() { };
	KStarHyperon( const vector< string >& args );

	enum UserVars { kPgamma = 0, kCosThetaX,  kCosThetaY,  kCosThetaZ, kPhi, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }

	string name() const { return "KStarHyperon"; }
    
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

  AmpParameter alpha;
  //AmpParameter Sigma;
  AmpParameter rho111;
  AmpParameter rho001; 
  AmpParameter Ox;
  AmpParameter P;
  AmpParameter T;
  AmpParameter Oz;
  AmpParameter polAngle;

  double polFraction;
  TH1D *polFrac_vs_E;

  std::vector<int> kIndices;
  std::vector<int> lambdaIndices;

};

#endif
