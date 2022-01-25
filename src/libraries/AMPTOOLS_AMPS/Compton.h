#if !defined(COMPTON)
#define COMPTON

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

class Compton : public UserAmplitude< Compton >
{
    
public:
	
	Compton() : UserAmplitude< Compton >() { };
	Compton( const vector< string >& args );
	
	string name() const { return "Compton"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:

	GDouble polAngle;
	GDouble polFraction;
    bool polInTree;
	TH1D *polFrac_vs_E;
};

#endif
