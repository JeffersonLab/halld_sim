#include "EtaPi0PlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

EtaPi0PlotGenerator::EtaPi0PlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  bookHistogram( khm12, new Histogram1D( 60, 0.0, 5.0, "hm12", "Mass( 1 2 )" ) );
  bookHistogram( khm13, new Histogram1D( 60, 0.0, 5.0, "hm13", "Mass( 1 3 )" ) );
  bookHistogram( khm23, new Histogram1D( 60, 0.0, 2.1, "hm23", "Mass( 2 3 )" ) );
  bookHistogram( kdltz, new Histogram2D( 80, 0.0, 25.0, 80, 0.0, 9.0, "dltz", "Dalitz Plot" ) );
  bookHistogram( cosT, new Histogram1D( 60, -1.1, 1.1, "cosT", "CosTheta") );
  bookHistogram( phiAng, new Histogram1D(40, -3.2, 3.2, "phiAng", "#phi") );
  bookHistogram( cosT_lab, new Histogram1D( 60, -1.1, 1.1, "cosT_lab", "CosThetaLab") );
  bookHistogram( phiAng_lab, new Histogram1D(40, -3.2, 3.2, "phiAng_lab", "#phi_{lab}") );
  bookHistogram( cosT_m23_lab, new Histogram2D(200, 0.6, 2.6, 100, -1.0, 1.0, "cosTLab_m23", "cos(#theta_{lab}) vs. Mass(#eta #pi^{0})" ) );
  bookHistogram( phi_m23_lab, new Histogram2D(200, 0.6, 2.6, 100, -3.2, 3.2, "PhiLab_m23", "#phi_{lab} vs. Mass(#eta #pi^{0})" ) );
  bookHistogram( PhiT, new Histogram1D(40, -3.2, 3.2, "PhiT", "#Phi") );
  bookHistogram( Omega, new Histogram1D(100, -3.2, 3.2, "Omega", "#Omega") );
  bookHistogram( cosT_m23, new Histogram2D(100, 1.8, 3.2, 100, -1.0, 1.0, "cosT_m23", "cos(#theta) vs. Mass(#eta #pi^{0})" ) );
  bookHistogram( cosT_phi, new Histogram2D(40, -3.2, 3.2, 100, -1.0, 1.0, "cosT_phi", "cos(#theta) vs. #phi" ) );
  bookHistogram( cosT_Phi, new Histogram2D(40, -3.2, 3.2, 100, -1.0, 1.0, "cosT_Phi", "cos(#theta) vs. #Phi" ) );

  polAngle = stod(results.configInfo()->amplitudeList("","","S0-").at(0)->factors().at(0).at(5));
}

void EtaPi0PlotGenerator::projectEvent( Kinematics* kin ){

  TLorentzVector P0 = kin->particle(0); //beam
  TLorentzVector P1 = kin->particle(1); //proton
  TLorentzVector P2 = kin->particle(2); //eta
  TLorentzVector P3 = kin->particle(3); //pi0

  fillHistogram( khm12, (P1+P2).M() );
  fillHistogram( khm13, (P1+P3).M() );
  fillHistogram( khm23, (P2+P3).M() );
  fillHistogram( kdltz, (P1+P2).M2(), (P2+P3).M2() );

  TLorentzVector resonance = P2 + P3;
  TLorentzRotation resRestBoost( -resonance.BoostVector() );

  TLorentzVector beam_res   = resRestBoost * P0;
  TLorentzVector recoil_res = resRestBoost * P1;
  TLorentzVector p2_res = resRestBoost * P2;

  // Angles of etapi system in LAB frame:
  GDouble locCosT_lab = resonance.CosTheta();
  GDouble locPhi_lab = resonance.Phi();

  fillHistogram( cosT_lab, locCosT_lab );
  fillHistogram( phiAng_lab, locPhi_lab );
  fillHistogram( cosT_m23_lab, resonance.M(), locCosT_lab );
  fillHistogram( phi_m23_lab, resonance.M(), locPhi_lab );
  
  // Helicity Frame:
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 y = (P0.Vect().Unit().Cross(-P1.Vect().Unit())).Unit();

  TVector3 x = y.Cross(z);

  TVector3 angles( (p2_res.Vect()).Dot(x),
                   (p2_res.Vect()).Dot(y),
                   (p2_res.Vect()).Dot(z) );

  Double_t cosTheta = angles.CosTheta();
  Double_t phi = angles.Phi();
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), P0.Vect().Unit().Dot(eps.Cross(y)));
  
  TVector3 eta = P2.Vect();
  GDouble omega = atan2(y.Dot(eta), P0.Vect().Unit().Dot(eta.Cross(y)));
  

  fillHistogram( PhiT, Phi); 
  fillHistogram( cosT, cosTheta);
  fillHistogram( phiAng, phi);
  fillHistogram( Omega, omega);
  fillHistogram( cosT_m23, (P2+P3).M(), cosTheta);
  fillHistogram( cosT_phi, phi, cosTheta);
  fillHistogram( cosT_Phi, Phi, cosTheta);
}
