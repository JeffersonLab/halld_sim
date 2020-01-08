//July 26th 2019, Based on DOI: 10.1016/0550-3213(84)90382-1
#if !defined(OMEGAPIANGAMP)
#define OMEGAPIANGAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"

#ifdef GPU_ACCELERATION
void GPUomegapiAngAmp_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
			    GDouble m_1p, GDouble m_w_1p, GDouble m_n_1p, GDouble m_1m, GDouble m_w_1m, GDouble m_n_1m,
                            GDouble m_0m, GDouble m_w_0m, GDouble m_n_0m, GDouble m_ds_ratio, GDouble m_phi0_1p,
                            GDouble m_theta_1p, GDouble m_phip_1p, GDouble m_phim_1p, GDouble m_psi_1p,
			    GDouble m_phi0_1m, GDouble m_theta_1m, GDouble m_phip_1m, GDouble m_phim_1m,
                            GDouble m_psi_1m, GDouble m_phi0_0m, GDouble m_theta_0m, bool useCutoff, GDouble polAngle, GDouble polFraction );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class omegapiAngAmp : public UserAmplitude< omegapiAngAmp >
{

public:
  
  omegapiAngAmp() : UserAmplitude< omegapiAngAmp >() { }
  omegapiAngAmp( const vector< string >& args );
  ~omegapiAngAmp(){}
  
  string name() const { return "omegapiAngAmp"; }
  
complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
  
  // **********************
  // The following lines are optional and can be used to precalcualte
  // user-defined data that the amplitudes depend on.
  
  // Use this for indexing a user-defined data array and notifying
  // the framework of the number of user-defined variables.
  
  enum UserVars { uv_Phi = 0, uv_Pgamma = 0, uv_mx = 0, uv_moment0 = 0
  , uv_moment1 = 0
  , uv_moment2 = 0
  , uv_moment3 = 0
  , uv_moment4 = 0
  , uv_moment5 = 0
  , uv_moment6 = 0
  , uv_moment7 = 0
  , uv_moment8 = 0
  , uv_moment9 = 0
  , uv_moment10 = 0
  , uv_moment11 = 0
  , uv_moment12 = 0
  , uv_moment13 = 0
  , uv_moment14 = 0
  , uv_moment15 = 0
  , uv_moment16 = 0
  , uv_moment17 = 0
  , uv_moment18 = 0
  , uv_moment19 = 0
  , uv_moment20 = 0
  , uv_moment21 = 0
  , uv_moment22 = 0
  , uv_moment23 = 0
  , uv_moment24 = 0
  , uv_calpha0 = 0
  , uv_calpha1 = 0
  , uv_calpha2 = 0
  , uv_calpha3 = 0
  , uv_calpha4 = 0
  , uv_calpha5 = 0
  , uv_calpha6 = 0
  , uv_calpha7 = 0
  , uv_calpha8 = 0
  , uv_calpha9 = 0
  , uv_calpha10 = 0
  , uv_calpha11 = 0
  , uv_calpha12 = 0
  , uv_calpha13 = 0
  , uv_calpha14 = 0
  , uv_calpha15 = 0
  , uv_calpha16 = 0
  , uv_calpha17 = 0
  , uv_calpha18 = 0
  , uv_calpha19 = 0
  , uv_calpha20 = 0
  , uv_calpha21 = 0
  , uv_calpha22 = 0
  , uv_calpha23 = 0
  , uv_calpha24 = 0,kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  
  // This function needs to be defined -- see comments and discussion
  // in the .cc file.
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
  
  // This is an optional addition if the calcAmplitude routine
  // can run with only the user-defined data and not the original
  // four-vectors.  It is used to optimize memory usage in GPU
  // based fits.
  bool needsUserVarsOnly() const { return true; }
  // **  end of optional lines **
  
 void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:
  bool useCutoff;
  AmpParameter m_1p;
  AmpParameter m_w_1p;
  AmpParameter m_n_1p;
  AmpParameter m_phi0_1p;
  AmpParameter m_phip_1p;
  AmpParameter m_phim_1p;
  AmpParameter m_theta_1p;
  AmpParameter m_psi_1p;

  AmpParameter m_1m;
  AmpParameter m_w_1m;
  AmpParameter m_n_1m;
  AmpParameter m_phi0_1m;
  AmpParameter m_phip_1m;
  AmpParameter m_phim_1m;
  AmpParameter m_theta_1m;
  AmpParameter m_psi_1m;

  AmpParameter m_0m;
  AmpParameter m_w_0m;
  AmpParameter m_n_0m;
  AmpParameter m_phi0_0m;
  AmpParameter m_theta_0m;

  AmpParameter m_ds_ratio;

  double polAngle, polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;

};

#endif
