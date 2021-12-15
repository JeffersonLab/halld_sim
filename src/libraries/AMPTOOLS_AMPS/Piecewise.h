#if !defined(PIECEWISE)
#define PIECEWISE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

//#ifdef GPU_ACCELERATION
//void GPUBreitWigner_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
//                          GDouble mass0, GDouble width0, int orbitL,
//                          int daught1, int daught2 );
//
//#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class Piecewise : public UserAmplitude< Piecewise >
{
  
public:
	
	Piecewise() : UserAmplitude< Piecewise >() {}
	Piecewise( const vector< string >& args );
	
  ~Piecewise(){}
  
	string name() const { return "Piecewise"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
  void init();
    
//#ifdef GPU_ACCELERATION
//
//  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
//
//	bool isGPUEnabled() const { return true; }
//
//#endif // GPU_ACCELERATION
  
private:
	
  float m_massMin, m_massMax;
  int m_daughter1, m_daughter2, m_nBins;
  complex<GDouble> one;
  complex<GDouble> zero;

  vector<AmpParameter> m_params1;
  vector<AmpParameter> m_params2;
  double m_width;
  string m_suffix;
  AmpParameter paramTest;
};

#endif
