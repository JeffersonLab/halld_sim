#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Polynomial.h"

Polynomial::Polynomial( const vector< string >& args ) :
UserAmplitude< Polynomial >( args )
{
  assert( args.size() == 8 );

        m_daughters = pair< string, string >( args[0], args[1] );

  m_mag_p0 = AmpParameter( args[2] );
  m_mag_p1 = AmpParameter( args[3] );
  m_mag_p2 = AmpParameter( args[4] );
  m_phase_p0 = AmpParameter( args[5] );
  m_phase_p1 = AmpParameter( args[6] );
  m_phase_p2 = AmpParameter( args[7] );

  registerParameter( m_mag_p0 );
  registerParameter( m_mag_p1 );
  registerParameter( m_mag_p2 );
  registerParameter( m_phase_p0 );
  registerParameter( m_phase_p1 );
  registerParameter( m_phase_p2 );
}

complex< GDouble >
Polynomial::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp;

  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){

    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P1 += Ptemp;
    Ptot += Ptemp;
  }

  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){

    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P2 += Ptemp;
    Ptot += Ptemp;
  }

  complex<GDouble> mass(Ptot.M(),0.0);

  double mag = m_mag_p0 + m_mag_p1*Ptot.M() + m_mag_p2*Ptot.M()*Ptot.M();
  double phase = m_phase_p0 + m_phase_p1*Ptot.M() + m_phase_p2*Ptot.M()*Ptot.M();

  return polar(mag,phase);
}
