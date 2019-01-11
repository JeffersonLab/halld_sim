#if !defined(TWOPIWT_SIGMA)
#define TWOPIWT_SIGMA

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

class TwoPiWt_sigma : public UserAmplitude< TwoPiWt_sigma >
{
  
public:
	
	TwoPiWt_sigma() : UserAmplitude< TwoPiWt_sigma >() {}
	TwoPiWt_sigma( const vector< string >& args );
	
	~TwoPiWt_sigma(){}
  
	string name() const { return "TwoPiWt_sigma"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_par1;    // for the moment assume W cross section has 2 parameters
  AmpParameter m_par2;
  
  pair< string, string > m_daughters;  
};

#endif
