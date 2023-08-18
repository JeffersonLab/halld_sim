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

#ifdef GPU_ACCELERATION
void
GPUSinglePS_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_r, int m_s );
#endif

class Kinematics;

class SinglePS : public UserAmplitude< SinglePS >
{
    
public:
	
	SinglePS() : UserAmplitude< SinglePS >() { };
	SinglePS( const vector< string >& args );
	SinglePS( int m_r, int m_s );
	
	string name() const { return "SinglePS"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

	// **********************
	// The following lines are optional and can be used to precalcualte
	// user-defined data that the amplitudes depend on.
	
	// Use this for indexing a user-defined data array and notifying
	// the framework of the number of user-defined variables.
	
	enum UserVars { uv_prod_Phi = 0, uv_beam_polFraction = 1, uv_beam_polAngle = 2, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }
	
	// This function needs to be defined -- see comments and discussion
	// in the .cc file.
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
	
	// This is an optional addition if the calcAmplitude routine
	// can run with only the user-defined data and not the original
	// four-vectors.  It is used to optimize memory usage in GPU
	// based fits.
	bool needsUserVarsOnly() const { return true; }
	
	// This is an optional addition if the UserVars are the same for each 
	// instance of an amplitude.  If it is not used, the memory footprint
	// grows dramatically as UserVars values are stored for each instance
	// of the amplitude.  NOTE: To use this make sure that UserVars only 
	// depend on kinematics and no arguments provided to the amplitude!
	bool areUserVarsStatic() const { return true; }

	void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

        bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
	
private:
        
	int m_r;
	int m_s;

	
	//AmpParameter polAngle;
	double polFraction;
	double polAngle;
	bool m_polInTree;
	TH1D *polFrac_vs_E;
};

#endif
