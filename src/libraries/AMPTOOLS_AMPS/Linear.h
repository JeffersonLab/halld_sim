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

  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;

  enum UserVars { uv_mass = 0, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }

  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

  bool needsUserVarsOnly() const { return false; }
  bool areUserVarsStatic() const { return false; }

  void updatePar( const AmpParameter& par );

private:

  pair< string, string > m_daughters;
  AmpParameter m_real_p0;
  AmpParameter m_real_p1;
  AmpParameter m_imag_p0;

  double imag_p1;

};

#endif
