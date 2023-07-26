#if !defined(LOWERVERTEXDELTA)
#define LOWERVERTEXDELTA

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
void GPULowerVertexDelta_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_d, int m_p, int m_c, int m_s );
#endif // GPU_ACCELERATION


class Kinematics;

class LowerVertexDelta : public UserAmplitude< LowerVertexDelta >
{
    
public:
	
	LowerVertexDelta() : UserAmplitude< LowerVertexDelta >() { };
	LowerVertexDelta( const vector< string >& args );
	LowerVertexDelta( int m_d, int m_p, int m_c, int m_s );

	enum UserVars { kCosTheta = 0, kPhi = 1, kNumUserVars };
	unsigned int numUserVars() const { return kNumUserVars; }
	
	string name() const { return "LowerVertexDelta"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
	void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

	bool needsUserVarsOnly() const { return true; }

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	int m_d;
	int m_p;
	int m_c;  
	int m_s;

	string lowerVertex;
};

#endif
