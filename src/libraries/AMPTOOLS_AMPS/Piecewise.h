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

#ifdef GPU_ACCELERATION
void GPUPiecewise_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* paramsRe, GDouble* paramsIm, int nBins, bool represReIm );
#endif // GPU_ACCELERATION

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
  
  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars  ) const;
	  
  // **********************
  // The following lines are optional and can be used to precalcualte
  // user-defined data that the amplitudes depend on.

  // Use this for indexing a user-defined data array and notifying
  // the framework of the number of user-defined variables.

  enum UserVars { uv_imassbin = 0, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }

  // This function needs to be defined -- see comments and discussion
  // in the .cc file.
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

  // This is an optional addition if the calcAmplitude routine
  // can run with only the user-defined data and not the original
  // four-vectors.  It is used to optimize memory usage in GPU
  // based fits.
  bool needsUserVarsOnly() const { return true; }

  // This is an optional addition if the UserVars are the same for each 
  // instance of an amplitude.  If it is not used, the memory footprint
  // grows dramatically as UserVars values are stored for each instance
  // of the amplitude.  NOTE: To use this make sure that UserVars only 
  // depend on kinematics and no arguments provided to the amplitude!
  bool areUserVarsStatic() const { return true; }

  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

  bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  float m_massMin, m_massMax;
  int m_nBins;
  string m_daughters;  
  complex<GDouble> one;
  complex<GDouble> zero;

  vector<AmpParameter> m_params1;
  vector<AmpParameter> m_params2;
  double m_width;
  string m_suffix;
  AmpParameter paramTest;
  bool m_represReIm;
};

#endif
