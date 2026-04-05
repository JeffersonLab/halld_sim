#if !defined(YLMREIM)
#define YLMREIM

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

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
// s=-1 indicates complex conjugation of Ylm

class Kinematics;

class Ylm_reim : public UserAmplitude< Ylm_reim >
{
    
public:
	
	Ylm_reim() : UserAmplitude< Ylm_reim >() { };
	Ylm_reim( const vector< string >& args );
	
	string name() const { return "Ylm_reim"; }
    
        complex<GDouble> calcAmplitude( GDouble** pKin ) const;
	
private:
        
  int m_j;
  int m_m;
  int m_r;// real or imaginary

  int m_phaseFactor;
	
};

#endif
