
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Ylm.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

Ylm::Ylm( const vector< string >& args ) :
UserAmplitude< Ylm >( args )
{
  assert( args.size() == 3 );
  
  m_j = atoi( args[0].c_str() );
  m_m = atoi( args[1].c_str() );
  m_s = atoi( args[2].c_str() );

  // make sure values are reasonable
  assert( abs( m_s ) == 1 );
  assert( abs( m_m ) <= m_j );

  m_phaseFactor = 1;
  // for complex conjugates: Ylm* = (-1)^m Yl-m
  if ( m_s < 0 ){
    m_phaseFactor = ( m_m % 2 == 0 ? 1 : -1 );
    m_m *= -1;
  }
  
}


complex< GDouble >
Ylm::calcAmplitude( GDouble** pKin ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  
  TLorentzVector resonance = p1 + p2;
  
  TLorentzRotation resRestBoost( -resonance.BoostVector() );
  
  TLorentzVector beam_res   = resRestBoost * beam;
  TLorentzVector recoil_res = resRestBoost * recoil;
  TLorentzVector p1_res = resRestBoost * p1;
  
  // Helicity frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  // or GJ frame?
  // TVector3 z = beam_res.Vect().Unit();
  TVector3 y = recoil_res.Vect().Cross(z).Unit();
  TVector3 x = y.Cross(z);
  
  TVector3 angles( (p1_res.Vect()).Dot(x),
                   (p1_res.Vect()).Dot(y),
                   (p1_res.Vect()).Dot(z) );
  
  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();

  return complex< GDouble >( static_cast< GDouble>( m_phaseFactor ) * Y( m_j, m_m, cosTheta, phi ) );
}

