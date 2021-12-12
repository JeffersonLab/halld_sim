

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"

PhaseOffset::PhaseOffset( const vector< string >& args ) :
UserAmplitude< PhaseOffset >( args )
{
  
  assert( args.size() == 1 );	
  m_phase = AmpParameter( args[0] );

  // need to register any free parameters so the framework knows about them
  registerParameter( m_phase );
}

complex< GDouble >
PhaseOffset::calcAmplitude( GDouble** pKin ) const
{
  complex <GDouble> a = polar(1.0, GDouble(m_phase));
  return a;
}
