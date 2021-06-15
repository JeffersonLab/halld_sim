#if !defined(ETAPB_TDIST)
#define ETAPB_TDIST

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

class EtaPb_tdist : public UserAmplitude< EtaPb_tdist >
{
  
public:
	
	EtaPb_tdist() : UserAmplitude< EtaPb_tdist >() {}
	EtaPb_tdist( const vector< string >& args );
	
	~EtaPb_tdist(){}
  
	string name() const { return "EtaPb_tdist"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

	//void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter ThetaSigma;    // for the moment assume W cross section has 5 parameters
  AmpParameter Phase;
  AmpParameter Bgen;
};

#endif
