#if !defined(DBLREGGE_FASTPI)
#define DBLREGGE_FASTPI

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPUDblRegge_FastPi_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble S0, GDouble b_pi, int charge );
#endif //GPU_ACCELERATION


using std::complex;
using namespace std;

class Kinematics;

class DblRegge_FastPi : public UserAmplitude< DblRegge_FastPi >
{

public:

        DblRegge_FastPi() : UserAmplitude< DblRegge_FastPi >() { };
        DblRegge_FastPi( const vector< string >& args );

        string name() const { return "DblRegge_FastPi"; }

        enum UserVars {u_s12=0,u_s23=1,u_t1=2,u_t2=3,u_s=4,u_u3=5,u_beamM2=6, u_p1M2=7, u_p2M2=8, u_recoilM2=9,u_up1=10, u_up2=11, kNumUserVars };
        
        unsigned int numUserVars() const {return kNumUserVars; }
	

	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
        void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
        double CHGM(double A, double B, double X) const;
        std::complex<double> cgamma(std::complex<double> z,int OPT) const;
        void updatePar( const AmpParameter& par );
        std::complex<double> V12(double alp1, double alp2, double eta) const;
        std::complex<double> DoubleRegge(int tau[2], double s, double sip[2], double alpp[2]) const;
        std::complex<double> ampEtaPi0(double par, int hel[3],  double inv[5], double mass2[4]) const;

        bool needsUserVarsOnly() const { return true; }
	bool areUserVarsStatic() const { return true; }
       // void init();

#ifdef GPU_ACCELERATION
	void launchGPUKernel ( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }
#endif // GPU_ACCELERATION

private:
  	int j;
        int fast;
	int charge; // 0 for neutral, 1 for charged
        AmpParameter b_pi;
        AmpParameter S0;
};

#endif
                     
