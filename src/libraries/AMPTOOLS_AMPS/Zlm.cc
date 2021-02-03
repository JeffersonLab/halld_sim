
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

#include "UTILITIES/BeamProperties.h"

Zlm::Zlm( const vector< string >& args ) :
UserAmplitude< Zlm >( args )
{
  assert( args.size() == 6 || args.size() == 7);
  
  m_j = atoi( args[0].c_str() );
  m_m = atoi( args[1].c_str() );
  m_r = atoi( args[2].c_str() );
  m_s = atoi( args[3].c_str() );

  polAngle = AmpParameter( args[4] );
  registerParameter( polAngle );
  
  polFraction = atof(args[5].c_str());
  
  // BeamProperties configuration file
  if (polFraction == 0){
    TString beamConfigFile = args[5].c_str();
    BeamProperties beamProp(beamConfigFile);
    polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
  }

  // make sure values are reasonable
  assert( abs( m_m ) <= m_j );
  // m_r = +1 for real
  // m_r = -1 for imag
  assert( abs( m_r ) == 1 );
  // m_s = +1 for 1 + Pgamma
  // m_s = -1 for 1 - Pgamma
  assert( abs( m_s ) == 1 );
  

  // Default reference frame: Helicity frame
  refFrame_sel = HELI;

  if(args.size()==7) {
     if(args[6] == "GJ")
        refFrame_sel = GJ;
     else if(args[6] == "Adair") 
        refFrame_sel = ADAIR;
   }
}


complex< GDouble >
Zlm::calcAmplitude( GDouble** pKin ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  
  TLorentzVector resonance = p1 + p2;
  
  TLorentzRotation resRestBoost( -resonance.BoostVector() );
  
  TLorentzVector beam_res   = resRestBoost * beam;
  TLorentzVector recoil_res = resRestBoost * recoil;
  TLorentzVector p1_res = resRestBoost * p1;
  

  // Switch reference frame:
  //
  // DEFAULT: Set z axis as required in helicity frame
  TVector3 z = -1. * recoil_res.Vect().Unit();

  // Switch frame if option was set:
  if(refFrame_sel == GJ) {
    z = beam_res.Vect().Unit();
  } else if(refFrame_sel == ADAIR) {
    z = beam.Vect().Unit();
  }



  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

  TVector3 x = y.Cross(z);
  
  TVector3 angles( (p1_res.Vect()).Dot(x),
                   (p1_res.Vect()).Dot(y),
                   (p1_res.Vect()).Dot(z) );
  
  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();

  GDouble Pgamma;  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
  
  if(polFraction > 0.) { // for fitting with constant polarization 
    Pgamma = polFraction;
  }
  else{
    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
      Pgamma = 0.;
    }
    else Pgamma = polFrac_vs_E->GetBinContent(bin);
  }
  
  GDouble Factor = sqrt(1 + m_s * Pgamma);
  GDouble zlm = 0;
  complex< GDouble > rotateY = polar(1., -1.*Phi);
  if (m_r == 1)
    zlm = real(Y( m_j, m_m, cosTheta, phi ) * rotateY);
  if (m_r == -1)
    zlm = imag(Y( m_j, m_m, cosTheta, phi ) * rotateY);

  return complex< GDouble >( static_cast< GDouble>( Factor ) * zlm );
}

