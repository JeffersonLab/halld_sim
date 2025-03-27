#if !defined(PHASESHIFT)
#define PHASESHIFT

#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>


using std::complex;
using namespace std;

class Kinematics;


class PhaseShift : public UserAmplitude< PhaseShift >
{
public:

    PhaseShift() : UserAmplitude< PhaseShift >() { }

    PhaseShift( const vector< string >& args );

    ~PhaseShift(){}

    string name() const { return "PhaseShift"; }

    complex< GDouble > calcAmplitude( GDouble** pKin ) const;

private:

    pair< string, string > m_daughters;
    AmpParameter m_p0;
    AmpParameter m_p1;

};

#endif
