
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Ylm_reim.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

Ylm_reim::Ylm_reim( const vector< string >& args ) :
UserAmplitude< Ylm_reim >( args )
{
  assert( args.size() == 3 );
  
  m_j = atoi( args[0].c_str() );
  m_m = atoi( args[1].c_str() );
  m_r = atoi( args[2].c_str() );// 1 for real, -1 for imaginary
  
  // make sure values are reasonable
  assert( abs( m_m ) <= m_j );
  assert( abs( m_r ) == 1 );

  // was some phase factor stuff here, but seems redundant when you complex conjugate in the fit automatically
  // fully removed m_s, not needed for polarization fraction or complex conjugation
}


complex <GDouble> Ylm_reim::calcAmplitude( GDouble** pKin ) const{
  
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

// normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

  TVector3 x = y.Cross(z);
  
  TVector3 angles( (p1_res.Vect()).Dot(x),
                   (p1_res.Vect()).Dot(y),
                   (p1_res.Vect()).Dot(z) );
  
  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();
  GDouble ylm_re = 0;
  GDouble ylm_im = 0;
  if(m_r==1){
    ylm_re = real(Y( m_j, m_m, cosTheta, phi ));
  } else{
    ylm_re = imag(Y( m_j, m_m, cosTheta, phi ));
  }
  return complex<GDouble> (ylm_re,ylm_im);
}

