#if !defined(TWOPIW_EMPTYTGT)
#define TWOPIW_EMPTYTGT

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

class TwoPiW_emptytgt : public UserAmplitude< TwoPiW_emptytgt >
{
  
public:
	
	TwoPiW_emptytgt() : UserAmplitude< TwoPiW_emptytgt >() {}
	TwoPiW_emptytgt( const vector< string >& args );
	
	~TwoPiW_emptytgt(){}
  
	string name() const { return "TwoPiW_emptytgt"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

	//void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_par1;    // for the moment assume W cross section has 2 parameters
  AmpParameter m_par2;
  
  pair< string, string > m_daughters;  
};

#endif
