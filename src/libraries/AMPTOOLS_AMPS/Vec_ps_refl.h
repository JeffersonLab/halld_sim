#if !defined(VEC_PS_REFL)
#define VEC_PS_REFL

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
// r=+/-1 indicates real/imaginary part of Vec_ps_refl
// s=+/-1 multiplies with sqrt(1+/- P_gamma)

class Kinematics;

class Vec_ps_refl : public UserAmplitude< Vec_ps_refl >
{
    
public:
	
	Vec_ps_refl() : UserAmplitude< Vec_ps_refl >() { };
	Vec_ps_refl( const vector< string >& args );
	
	string name() const { return "Vec_ps_refl"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

	// **********************
	// The following lines are optional and can be used to precalcualte
	// user-defined data that the amplitudes depend on.
	
	// Use this for indexing a user-defined data array and notifying
	// the framework of the number of user-defined variables.
	
	enum UserVars { uv_cosTheta = 0, uv_Phi = 1, uv_cosThetaH = 2, uv_PhiH = 3, uv_prod_angle = 4, uv_Pgamma = 5, uv_dalitz_z = 6, uv_dalitz_sin3theta = 7, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }
	
	// This function needs to be defined -- see comments and discussion
	// in the .cc file.
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
	
	// This is an optional addition if the calcAmplitude routine
	// can run with only the user-defined data and not the original
	// four-vectors.  It is used to optimize memory usage in GPU
	// based fits.
	bool needsUserVarsOnly() const { return true; }

	void updatePar( const AmpParameter& par );
	
private:
        
	int m_j;
	int m_m;
	int m_l;
	int m_r;
	int m_s;
	int m_3pi;
	
	AmpParameter dalitz_alpha;
	AmpParameter dalitz_beta;
	AmpParameter dalitz_gamma;
	AmpParameter dalitz_delta;
	
	AmpParameter polAngle;
	
	double polFraction;
	TH1D *polFrac_vs_E;
};

#endif
