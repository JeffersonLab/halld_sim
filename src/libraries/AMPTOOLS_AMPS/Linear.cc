#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Linear.h"

Linear::Linear( const vector< string >& args ) :
UserAmplitude< Linear >( args )
{
  assert( args.size() == 5 );

  m_daughters = pair< string, string >( args[0], args[1] );

  m_real_p0 = AmpParameter( args[2] );
  m_real_p1 = AmpParameter( args[3] );
  m_imag_p0 = AmpParameter( args[4] );

  registerParameter( m_real_p0 );
  registerParameter( m_real_p1 );
  registerParameter( m_imag_p0 );
}

complex< GDouble >
Linear::calcAmplitude( GDouble** pKin ) const
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

  double imag_p1 = m_real_p1 * m_imag_p0 / m_real_p0;

  double real_tot = m_real_p0 + m_real_p1*Ptot.M();
  double imag_tot = m_imag_p0 + imag_p1*Ptot.M();

  complex< GDouble > ans( real_tot, imag_tot );
  return ans;
}
