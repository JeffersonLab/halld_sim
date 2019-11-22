#if !defined(TWOPITDIST)
#define TWOPITDIST

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

class TwoPitdist : public UserAmplitude< TwoPitdist >
{
  
public:
	
	TwoPitdist() : UserAmplitude< TwoPitdist >() {}
	TwoPitdist( const vector< string >& args );
	
	~TwoPitdist(){}
  
	string name() const { return "TwoPitdist"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter Bslope;    // for the moment assume W cross section has 4 parameters
  AmpParameter Bgen;
  
  pair< string, string > m_daughters;  
};

#endif
