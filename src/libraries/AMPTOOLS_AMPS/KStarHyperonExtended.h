#if !defined(KSTARHYPERONEXTENDED)
#define KSTARHYPERONEXTENDED

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include <string>
#include <complex>
#include <bitset>
#include <vector>

using std::complex;
using namespace std;

struct term_desc {
	std::string name;   
	int i;
	int j;
	int k;
};

struct term_ijk {
	int i;
	int j;
	int k;
};
  
#ifdef GPU_ACCELERATION
void GPUKStarHyperonExtended_exec(
    dim3 dimGrid,
    dim3 dimBlock,
    GPU_AMP_PROTO,
    int model,
    GDouble alpha,
    GDouble* coeffs,                 // host coeff array (size = (model==2 ? 39 : nTerms))  
    term_ijk* terms,                 // host (i,j,k) array (size nTerms; nullptr if nTerms==0) 
    int nTerms,                      // number of (i,j,k) terms; pass 0 for model==2        
    bool schillingIncluded,
    GDouble* rhos,                   // size 9 if schillingIncluded else nullptr
    GDouble polAngle
);
#endif



class Kinematics;

class KStarHyperonExtended : public UserAmplitude< KStarHyperonExtended >
{
public:

    KStarHyperonExtended() : UserAmplitude< KStarHyperonExtended >() { };
    KStarHyperonExtended( const vector< string >& args );

    enum UserVars {
        // g_func cache (0–1)
        kPgamma = 0,
        kBigPhi,

        // n_func cache (2–4)
        kN0,                                 
        kN1,                                 
        kN2,                                 

        // b_func cache (5–10)
        kB0,                                 
        kB1,                                 
        kB2,                                 
        kB3,                                 
        kB4,                                 
        kB5,                                 

        kNumUserVars
    };
    unsigned int numUserVars() const { return kNumUserVars; }

    string name() const { return "KStarHyperonExtended"; }

    complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
    void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

    // we can calculate everything we need from userVars so allow the framework to purge the four-vectors
    bool needsUserVarsOnly() const { return true; }

    // user vars depend only on event kinematics (and beam energy for pol fraction), not on fit params
    bool areUserVarsStatic() const { return true; }

#ifdef GPU_ACCELERATION
    void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
    bool isGPUEnabled() const { return true; }
#endif // GPU_ACCELERATION

private:

    AmpParameter alpha;
	std::vector<AmpParameter> coeffs;             
	std::vector<AmpParameter> rhos; 
	AmpParameter polAngle;
          
	std::bitset<54> termMask;   // 54-bit mask encoding which (i,j,k) terms are active
    std::vector<int> iList;                 
    std::vector<int> jList;                 
    std::vector<int> kList;                 

    double polFraction;
    TH1D* polFrac_vs_E;

    std::vector<int> kIndices;
    std::vector<int> lambdaIndices;

    int model;
	bool schillingIncluded;     

};

#endif
