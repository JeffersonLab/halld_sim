#if !defined(PI0REGGE)
#define PI0REGGE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class PiPlusRegge : public UserAmplitude< PiPlusRegge >
{
    
public:
	
	PiPlusRegge() : UserAmplitude< PiPlusRegge >() { };
	PiPlusRegge( const vector< string >& args );
	
	string name() const { return "PiPlusRegge"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:

	GDouble PolPlane;

	TH1D *totalFlux_vs_E;
	TH1D *polFlux_vs_E;
	TH1D *polFrac_vs_E;
};

#endif
