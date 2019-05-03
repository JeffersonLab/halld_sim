#if !defined(LAMBDA1520TDIST)
#define LAMBDA1520TDIST

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

class Lambda1520tdist : public UserAmplitude< Lambda1520tdist >
{
  
public:
	
	Lambda1520tdist() : UserAmplitude< Lambda1520tdist >() {}
	Lambda1520tdist( const vector< string >& args );
	
	~Lambda1520tdist(){}
  
	string name() const { return "Lambda1520tdist"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter Bslope;    // for the moment assume W cross section has 4 parameters
  AmpParameter exponent;
  AmpParameter Bgen;
  
  pair< string, string > m_daughters;  
};

#endif
