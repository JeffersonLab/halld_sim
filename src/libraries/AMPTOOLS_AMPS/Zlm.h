#if !defined(ZLM)
#define ZLM

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

// A class for describing the angular portion of the decay R->12
// in the reflectivity basis for photon beams: m can be negative!
// particles 1 and 2 are pseudoscalars
//
// j,m are the total and z projection of the spin of R
// r=+/-1 indicates real/imaginary part of Zlm
// s=+/-1 multiplies with sqrt(1+/- P_gamma)

class Kinematics;

class Zlm : public UserAmplitude< Zlm >
{
    
public:
	
	Zlm() : UserAmplitude< Zlm >() { };
	Zlm( const vector< string >& args );
	
	string name() const { return "Zlm"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:
        
  int m_j;
  int m_m;
  int m_r;
  int m_s;
	
  AmpParameter polAngle;

  double polFraction;
  TH1D *polFrac_vs_E;
};

#endif
