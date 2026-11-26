#if !defined(INVGAUSS)
#define INVGAUSS

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;  

class Kinematics;

class InvGauss : public UserAmplitude< InvGauss >
{
  
public:
	
  InvGauss() : UserAmplitude< InvGauss >() {}
  InvGauss( const vector< string >& args );	
  ~InvGauss(){}

  string name() const { return "InvGauss"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;	  

  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

	//void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
  
  AmpParameter Lambda;
  AmpParameter Mu;
  
  string m_daughters;  
};

#endif
