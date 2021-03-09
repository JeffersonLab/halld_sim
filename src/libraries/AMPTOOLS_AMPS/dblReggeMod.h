#if !defined(DBLREGGEMOD)
#define DBLREGGEMOD

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class dblReggeMod : public UserAmplitude< dblReggeMod >
{

public:

	dblReggeMod() : UserAmplitude< dblReggeMod >() { };
        dblReggeMod( const vector< string >& args );

        string name() const { return "dblReggeMod"; }

        complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	double CHGM(double A, double B, double X) const;
	std::complex<double> cgamma(std::complex<double> z,int OPT) const;
	void updatePar( const AmpParameter& par );
	std::complex<double> V12(double alp1, double alp2, double eta) const;
	std::complex<double> DoubleRegge(int tau[2], double s, double si[2], double alp[2]) const;
	std::complex<double> ampEtaPi0(double par[4], int hel[3],  double inv[5], double mass2[4]) const;


private:

//	double polFraction;
        int j;	
	int fast;
	AmpParameter b_eta, b_pi, a_eta, a_pi;
	AmpParameter S0;
};

#endif
