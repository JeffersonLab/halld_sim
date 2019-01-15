//april 9th 2018, chung 1975
#if !defined(OMEGAPIANGAMP)
#define OMEGAPIANGAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"

#include "GPUManager/GPUCustomTypes.h"

using std::complex;
using namespace std;


#ifdef GPU_ACCELERATION
void
GPUomegapiAngAmp_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		   int polBeam, GDouble polFrac,
		   int J_X, int Par_X, int L_X, int I_X, int epsilon_R, 
		   int Iz_b1, int Iz_pi,
		   GDouble u_rho_1, GDouble u_rho_3, 
		   GDouble u_omega_1, GDouble u_omega_3,
		   GDouble u_b1_0, GDouble u_b1_2, 
		   GDouble G0_omega, GDouble G0_b1, bool orthocheck);
#endif



class Kinematics;

class omegapiAngAmp : public UserAmplitude< omegapiAngAmp >
{

public:
  
  omegapiAngAmp() : UserAmplitude< omegapiAngAmp >() { }
  omegapiAngAmp( const vector< string >& args );
  
  string name() const { return "omegapiAngAmp"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

  GDouble u_rho(int J_rho) const;
  GDouble u_omega(int L_omega) const;
  GDouble u_b1(int L_b1) const;

  inline complex<GDouble> BreitWigner(GDouble m0, GDouble Gamma0, int L,
			       TLorentzVector &P1,TLorentzVector &P2) const;
  inline GDouble CB(int j1, int j2, int m1, int m2, int J, int M) const;

 
private:
  
  int mpolBeam;
  AmpParameter mpolFrac;
  int mJ_X, mPar_X, mL_X, mI_X, mepsilon_R;
  
  GDouble m_u_rho_1, m_u_rho_3;
  GDouble m_u_omega_1, m_u_omega_3;
  GDouble m_u_b1_0, m_u_b1_2, mG0_omega, mG0_b1;
  bool m_ORTHOCHECK, m_fastCalc;
  bool m_disableBW_omega, m_disableBW_b1;

  vector< int > mIz;
  

  AmpParameter polAngle;

  GDouble m_1p;
  GDouble w_1p;
  GDouble n_1p;
  GDouble phi0_1p;
  GDouble phip_1p;
  GDouble phim_1p;
  GDouble theta_1p;
  GDouble psi_1p;
  
  GDouble m_1m;
  GDouble w_1m;
  GDouble n_1m;
  GDouble phi0_1m;
  GDouble phip_1m;
  GDouble phim_1m;
  GDouble theta_1m;
  GDouble psi_1m;

  GDouble m_0m;
  GDouble w_0m;
  GDouble n_0m;
  GDouble phi0_0m;
  GDouble theta_0m;

  GDouble ds_ratio;

  double polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;



#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
  bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  

};

#endif
