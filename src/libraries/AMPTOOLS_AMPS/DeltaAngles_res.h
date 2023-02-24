#if !defined(DELTAANGLES_RES)
#define DELTAANGLES_RES

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

class DeltaAngles_res : public UserAmplitude< DeltaAngles_res >
{
    
public:
	
	DeltaAngles_res() : UserAmplitude< DeltaAngles_res >() { };
	DeltaAngles_res( const vector< string >& args );
	
	string name() const { return "DeltaAngles_res"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;

  
private:
  
	AmpParameter rho011;
	AmpParameter rho031;
	AmpParameter rho03m1;
	
	AmpParameter rho111;
	AmpParameter rho133;
	AmpParameter rho131;
	AmpParameter rho13m1;
	
	AmpParameter rho231;
	AmpParameter rho23m1;

	GDouble polFraction=0.;
	GDouble polAngle=-1;
	TH1D *polFrac_vs_E;

	string lowerVertex;
	string upperVertex;

};

#endif
