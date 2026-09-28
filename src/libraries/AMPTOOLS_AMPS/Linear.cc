#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Linear.h"

// this amplitude returns a function that is linear in
// invariant mass and has a constant phase.  It takes
// the form of e^(i*phi) ( a + b*M ), where M is an
// invaraiant mass and phi, a, and b, are all real
// numbers -- the computation below is equivalent but
// done in terms of real and imaginary parts

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
Linear::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  GDouble mass = userVars[kMass];

  double real_tot = m_real_p0 + m_real_p1 * mass;
  double imag_tot = m_imag_p0 + m_imag_p1 * mass;

  complex< GDouble > ans( real_tot, imag_tot );
  return ans;
}

void Linear::calcUserVars( GDouble** pKin, GDouble* userVars ) const
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
  userVars[kMass] = Ptot.M();
}

void
Linear::updatePar( const AmpParameter& par )
{
  m_imag_p1 = m_real_p1 * m_imag_p0 / m_real_p0;
}


#ifdef GPU_ACCELERATION
void
Linear::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  GPULinear_exec( dimGrid,  dimBlock, GPU_AMP_ARGS,
                  m_real_p0, m_real_p1, m_imag_p0, m_imag_p1 );

}
#endif //GPU_ACCELERATION

