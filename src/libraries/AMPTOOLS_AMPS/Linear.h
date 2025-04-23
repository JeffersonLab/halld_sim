#if !defined(LINEAR)
#define LINEAR

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


class Linear : public UserAmplitude< Linear >
{
public:

  Linear() : UserAmplitude< Linear >() { }

  Linear( const vector< string >& args );

  ~Linear(){}

  string name() const { return "Linear"; }

  complex< GDouble > calcAmplitude( GDouble** pKin ) const;


private:

  pair< string, string > m_daughters;
  AmpParameter m_real_p0;
  AmpParameter m_real_p1;
  AmpParameter m_imag_p0;

};

#endif
