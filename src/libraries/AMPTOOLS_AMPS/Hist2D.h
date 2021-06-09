#if !defined(HIST2D)
#define HIST2D

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class Hist2D : public UserAmplitude< Hist2D >
{
    
public:
	
	Hist2D() : UserAmplitude< Hist2D >() { };
	Hist2D( const vector< string >& args );
	
	string name() const { return "Hist2D"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:
	
        string fileName, histName, histType, particleList;
	TH2 *hist2D;
};

#endif
