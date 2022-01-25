#if !defined(LAMBDA1520ANGLES)
#define LAMBDA1520ANGLES

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

class Lambda1520Angles : public UserAmplitude< Lambda1520Angles >
{
    
public:
	
	Lambda1520Angles() : UserAmplitude< Lambda1520Angles >() { };
	Lambda1520Angles( const vector< string >& args );
	
	string name() const { return "Lambda1520Angles"; }
    
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
    bool polInTree;

};

#endif
