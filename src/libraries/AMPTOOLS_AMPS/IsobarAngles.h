#if !defined(ISOBARANGLES)
#define ISOBARANGLES

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class IsobarAngles : public UserAmplitude< IsobarAngles >
{

public:
	
	IsobarAngles() : UserAmplitude< IsobarAngles >() { }
	IsobarAngles( const vector< string >& args );
	
	string name() const { return "IsobarAngles"; }
	
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
private:
	
	int m_jX;
	int m_lX;	

	string m_daughtX;
	vector<string> m_daughtI;
	vector<int> m_jI;
	int m_nIsobars;
    
};

#endif
