#if !defined(TWOPIETAS_TDIST)
#define TWOPIETAS_TDIST

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

class TwoPiEtas_tdist : public UserAmplitude< TwoPiEtas_tdist >
{
  
public:
	
	TwoPiEtas_tdist() : UserAmplitude< TwoPiEtas_tdist >() {}
	TwoPiEtas_tdist( const vector< string >& args );
	
	~TwoPiEtas_tdist(){}
  
	string name() const { return "TwoPiEtas_tdist"; }
  
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
