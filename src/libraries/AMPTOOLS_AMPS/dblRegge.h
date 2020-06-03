#if !defined(DBLREGGE)
#define DBLREGGE

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

class dblRegge : public UserAmplitude< dblRegge >
{

public:

        dblRegge() : UserAmplitude< dblRegge >() { };
        dblRegge( const vector< string >& args );

        string name() const { return "dblRegge"; }

        complex< GDouble > calcAmplitude( GDouble** pKin ) const;

private:

        double polFraction;
        int j;
        int fast;
        double b;
};

#endif
