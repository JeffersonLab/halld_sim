#if !defined(TWOPINC_TDIST)
#define TWOPINC_TDIST

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

class TwoPiNC_tdist : public UserAmplitude< TwoPiNC_tdist >
{
  
public:
	
	TwoPiNC_tdist() : UserAmplitude< TwoPiNC_tdist >() {}
	TwoPiNC_tdist( const vector< string >& args );
	
	~TwoPiNC_tdist(){}
  
	string name() const { return "TwoPiNC_tdist"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter ThetaSigma;    // for the moment assume W cross section has 5 parameters
  AmpParameter Phase;
  AmpParameter Bgen;
  
  pair< string, string > m_daughters;  
};

#endif
