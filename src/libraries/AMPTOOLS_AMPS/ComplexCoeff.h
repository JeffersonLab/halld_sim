#if !defined(COMPLEXCOEFF)
#define COMPLEXCOEFF

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

class ComplexCoeff : public UserAmplitude< ComplexCoeff >
{  
public:
  
  ComplexCoeff() : UserAmplitude< ComplexCoeff >() { }
  ComplexCoeff( const vector< string >& args );
  
  ~ComplexCoeff(){}
  
  string name() const { return "ComplexCoeff"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

  void updatePar( const AmpParameter& par );

private:
	
  AmpParameter m_param1;
  AmpParameter m_param2;
  bool m_represReIm;

  complex< GDouble > m_value;
};

#endif
