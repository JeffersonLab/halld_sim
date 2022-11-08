

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"

ComplexCoeff::ComplexCoeff( const vector< string >& args ) :
UserAmplitude< ComplexCoeff >( args )
{
  
  assert( args.size() == 3 );	
  m_param1 = AmpParameter( args[0] );
  m_param2 = AmpParameter( args[1] );

  // switch between representation of complex parameters in Re/Im format and Mag/Phi format
  if(args[2] == "ReIm")
    m_represReIm = true;
  else if(args[2] == "MagPhi")
    m_represReIm = false;
  else
     cout << "ERROR: '" << args[2] << "' is not a defined mode for the Piecewise amplitude! Choose 'ReIm' or 'MagPhi'." << endl;

  // need to register any free parameters so the framework knows about them
  registerParameter( m_param1 );
  registerParameter( m_param2 );
}

complex< GDouble >
ComplexCoeff::calcAmplitude( GDouble** pKin ) const
{
  complex<GDouble> a(m_param1,m_param2);
  if(!m_represReIm)
    a = polar(fabs(GDouble(m_param1)), GDouble(m_param2));
  
  return a;
}
