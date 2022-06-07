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
	GDouble CHGM(GDouble A, GDouble B, GDouble X) const;
	std::complex<GDouble> cgamma(std::complex<GDouble> z,int OPT) const;
	void updatePar( const AmpParameter& par );
	std::complex<GDouble> V12(GDouble alp1, GDouble alp2, GDouble eta) const;
	std::complex<GDouble> DoubleRegge(int tau[2], GDouble s, GDouble si[2], GDouble alp[2]) const;
	std::complex<GDouble> ampEtaPi0(GDouble par[4], int hel[3],  GDouble inv[5], GDouble mass2[4]) const;


private:

//	GDouble polFraction;
        int j;	
	int fast;
	AmpParameter b_eta, b_pi, a_eta, a_pi;
	AmpParameter S0;
};

#endif
