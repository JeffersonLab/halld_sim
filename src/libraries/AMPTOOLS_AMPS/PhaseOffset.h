#if !defined(PHASEOFFSET)
#define PHASEOFFSET

#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>


using std::complex;
using namespace std;

class Kinematics;

#ifdef GPU_ACCELERATION
void GPUPhaseOffset_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO);

#endif


class PhaseOffset : public UserAmplitude< PhaseOffset >
{  
public:
  
  PhaseOffset() : UserAmplitude< PhaseOffset >() { }
  PhaseOffset( const vector< string >& args );
  
  ~PhaseOffset(){}
  
  string name() const { return "PhaseOffset"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
      
#ifdef GPU_ACCELERATION
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const{
    GPUPhaseOffset_exec(dimGrid, dimBlock, GPU_AMP_ARGS);
  };
  
  bool isGPUEnabled() const { return true; }
#endif

private:
	
  AmpParameter m_phase;

};

#endif
