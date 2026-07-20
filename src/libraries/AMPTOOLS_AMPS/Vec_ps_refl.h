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

#ifdef GPU_ACCELERATION
void
GPUVec_ps_refl_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO );
#endif

class Kinematics;

class Vec_ps_refl : public UserAmplitude< Vec_ps_refl >
{
    
public:
	
	Vec_ps_refl() : UserAmplitude< Vec_ps_refl >() { };
	Vec_ps_refl( const vector< string >& args );
	Vec_ps_refl( int m_j, int m_m, int m_l, int m_r, int m_s );
	
	string name() const { return "Vec_ps_refl"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

	// **********************
	// The following lines are optional and can be used to precalculate
	// user-defined data that the amplitudes depend on.
	
	// Use this for indexing a user-defined data array and notifying
	// the framework of the number of user-defined variables.
		
	//enum UserVars { uv_cosTheta = 0, uv_Phi, uv_cosThetaH, uv_PhiH,
	//                uv_prod_Phi, uv_MX, uv_MVec, uv_MPs, uv_beam_polFraction,
	//                uv_beam_polAngle, kNumUserVars };
	enum UserVars { uv_ampRe = 0, uv_ampIm, kNumUserVars };
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
	bool areUserVarsStatic() const { return false; }

	void updatePar( const AmpParameter& par );

	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//....oooOO0OOooo........ Structs ........oooOO0OOooo.....
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	struct VecPsReflArgs {
	  // Required positional args
	  int j;      // Total spin
	  int m;      // Spin projection
	  int l;      // Relative orbital angular momentum
	  int r;      // +1/-1 for (r)eal/imaginary part of the amplitude
	  int s;      // +1/-1 (s)ign of the P_gamma term
	
	  // Polarization information
	  bool   polInfoInPhotonP4    = true;
	  double polAngle             = -2.0;   // -2 = not set  
	  double polFraction          = -2.0;   // -2 = not set
	  // Amorphous sometimes is referred as having a polAngle -1, so we use -2 to 
	  // avoid possible confusion for both polAngle and polFraction
	  std::string polFile      = ""; // path to .root file with polarization vs E_beam
	  std::string polHist      = ""; // name of TH1D inside polFile
	
	  // Default is 2-body vector decay, but if omega is used as vector:
	  // Omega decay mode
	  bool omega3pi  = false;
	  bool omegagpi0 = false;
	  int  gHelicity = 0;  // +1/-1 helicity of the bachelor photon in omega->g pi0
	
	  bool noBarrier = false;    // Turn off barrier factors
	};
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//....oooOO0OOooo........ Helper Functions ........oooOO0OOooo.....
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	static bool isValidNumber(const std::string& argInput, double &value){
		char* end = nullptr;
		errno = 0;  // reset global error
		value = std::strtod(argInput.c_str(), &end);
	
		// Check if 
		// (1) no conversion was performed 
		// (2) there are leftover characters
		// (3) an overflow/underflow occurred   
		if(end == argInput.c_str() || *end != '\0' || errno != 0) {
			return false;  // not a valid number
		}
		// If "end" points to the end of string, it's fully numeric
		return true;
	}
	
	static double parseValidatedNumber(const std::string& label, 
									   const std::string& argInput,
									   const std::string& context = "Vec_ps_refl"){
		double tmpValue = 0.0;
		if(!isValidNumber(argInput, tmpValue)){
		  throw std::invalid_argument("[ " + context + " ] invalid " +
			                           label + ": " + argInput);
		}
		return tmpValue;
	}

	static VecPsReflArgs parsedArgs(const std::vector<std::string>& args, 
	  							    const std::string& context = "Vec_ps_refl");

#ifdef GPU_ACCELERATION

	void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

        bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
	
private:
        
	int m_j;
	int m_m;
	int m_l;
	int m_r;
	int m_s;
	
	bool m_polInfoInPhotonP4;
	GDouble m_polFraction;
	GDouble m_polAngle;
	TH1D *m_polFracVsE;

	bool m_3pi;
	bool m_gpi0;
	int  m_ghel;
	bool m_noBarrier;
	
	
};

#endif
