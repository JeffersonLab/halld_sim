#if !defined(TWOPIW_BROKENETAS)
#define TWOPIW_BROKENETAS

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

class TwoPiW_brokenetas : public UserAmplitude< TwoPiW_brokenetas >
{
  
public:
	
	TwoPiW_brokenetas() : UserAmplitude< TwoPiW_brokenetas >() {}
	TwoPiW_brokenetas( const vector< string >& args );
	
	~TwoPiW_brokenetas(){}
  
	string name() const { return "TwoPiW_brokenetas"; }
  
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
