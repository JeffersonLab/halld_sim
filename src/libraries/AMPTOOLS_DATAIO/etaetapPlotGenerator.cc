#include "AMPTOOLS_DATAIO/etaetapPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Ylm.h"


// particle 1 is proton, particle 2 is eta', particle 3 is eta
etaetapPlotGenerator::etaetapPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Ylm());

  bookHistogram( khm12, new Histogram1D( 50, 1.0, 4.2, "hm12", "Mass(p#eta') [GeV/c^{2}]" ) );
  bookHistogram( khm13, new Histogram1D( 50, 1.0, 4.0, "hm13", "Mass(p#eta) [GeV/c^{2}]" ) );
  bookHistogram( khm23, new Histogram1D( 20, 1.4, 3.4, "hm23", "M(#eta#eta') [GeV/c^{2}]" ) );
  bookHistogram( khm1, new Histogram1D( 100, 0.5, 1.5, "hm1", "Mass( p )" ) );
  bookHistogram( khm2, new Histogram1D( 100, 0.7, 1.5, "hm2", "Mass( #eta' )" ) );
  bookHistogram( khm3, new Histogram1D( 100, 0.0, 1.0, "hm3", "Mass( #eta )" ) );
  bookHistogram( kdltz, new Histogram2D( 80, 0.0, 25.0, 80, 0.0, 9.0, "dltz", "Dalitz Plot" ) );
  bookHistogram( cosT, new Histogram1D( 22, -1.1, 1.1, "cosT", "CosTheta") );
  bookHistogram( phiAng, new Histogram1D(25, -3.2, 3.2, "phiAng", "#phi") );
  bookHistogram( cosT_lab, new Histogram1D( 55, -1.1, 1.1, "cosT_lab", "CosThetaLab") );
  bookHistogram( phiAng_lab, new Histogram1D(50, -3.2, 3.2, "phiAng_lab", "#phi_{lab}") );
  bookHistogram( cosT_m23_lab, new Histogram2D(50, 0.6, 2.6, 50, -1.0, 1.0, "cosTLab_m23", "cos(#theta_{lab}) vs. Mass(#eta #eta')" ) );
  bookHistogram( phi_m23_lab, new Histogram2D(50, 0.6, 2.6, 50, -3.2, 3.2, "PhiLab_m23", "#phi_{lab} vs. Mass(#eta #eta')" ) );
  bookHistogram( PhiT, new Histogram1D(25, -3.2, 3.2, "PhiT", "#Phi") );
  bookHistogram( Omega, new Histogram1D(50, -3.2, 3.2, "Omega", "#Omega") );
  bookHistogram( cosT_m23, new Histogram2D(50, 0.7, 2.7, 50, -1.0, 1.0, "cosT_m23", "cos(#theta) vs. Mass(#eta #eta')" ) );
  bookHistogram( cosT_phi, new Histogram2D(40, -3.2, 3.2, 50, -1.0, 1.0, "cosT_phi", "cos(#theta) vs. #phi" ) );
  bookHistogram( cosT_Phi, new Histogram2D(40, -3.2, 3.2, 50, -1.0, 1.0, "cosT_Phi", "cos(#theta) vs. #Phi" ) );
}

void etaetapPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void etaetapPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  // obtain the polarzation angle for this event by getting the list of amplitudes
  // associated with this reaction -- we assume here that the first amplitude in the list is a Zlm amplitude
  // take the sixth argument of the first factor of the first amplitude in the first sum

  int nargs  = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).size();
  TVector3 eps;
  TLorentzVector P0;


  // was originally just for 7, but I think it works for 9 as well
  if(nargs==7||nargs==9) {
     double polAngle = stod(cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).at(6));
     polAngle = 0.0;
     eps.SetXYZ(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
     P0 = kin->particle(0);
  } else if(nargs==5) {
     P0.SetPxPyPzE(0., 0., kin->particle(0).E(),kin->particle(0).E());
     eps.SetXYZ(kin->particle(0).Px(),kin->particle(0).Py(), 0.);
  } else {
     std::cout << "Zlm does not have required number of arguments (4 OR 6, but found " << nargs-1 << ")! Aborting." << std::endl;
     std::cout << cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).at(0) << std::endl;
     return;
  }

  TLorentzVector P1 = kin->particle(1); //proton
  TLorentzVector P2 = kin->particle(2); //eta
  TLorentzVector P3 = kin->particle(3); //pi0

  fillHistogram( khm12, (P1+P2).M() );
  fillHistogram( khm13, (P1+P3).M() );
  fillHistogram( khm23, (P2+P3).M() );
  fillHistogram( khm1, (P1).M() );
  fillHistogram( khm2, (P2).M() );
  fillHistogram( khm3, (P3).M() );
  fillHistogram( kdltz, (P1+P2).M2(), (P2+P3).M2() );

  TLorentzVector resonance = P2 + P3;
  TLorentzRotation resRestBoost( -resonance.BoostVector() );

  TLorentzVector beam_res   = resRestBoost * P0;
  TLorentzVector recoil_res = resRestBoost * P1;
  TLorentzVector p3_res = resRestBoost * P3;

  // Angles of etapi system in LAB frame:
  GDouble locCosT_lab = resonance.CosTheta();
  GDouble locPhi_lab = resonance.Phi();

  fillHistogram( cosT_lab, locCosT_lab );
  fillHistogram( phiAng_lab, locPhi_lab );
  fillHistogram( cosT_m23_lab, resonance.M(), locCosT_lab );
  fillHistogram( phi_m23_lab, resonance.M(), locPhi_lab );

  // Helicity Frame:
  TVector3 z = -1. * recoil_res.Vect().Unit();
  // or GJ frame?
//   TVector3 z = beam_res.Vect().Unit();
  
   TVector3 y = (P0.Vect().Unit().Cross(-P1.Vect().Unit())).Unit();

  TVector3 x = y.Cross(z);

  TVector3 angles( (p3_res.Vect()).Dot(x),
                   (p3_res.Vect()).Dot(y),
                   (p3_res.Vect()).Dot(z) );

  Double_t cosTheta = angles.CosTheta();
  Double_t phi = angles.Phi();
  
  GDouble Phi = atan2(y.Dot(eps), P0.Vect().Unit().Dot(eps.Cross(y)));
  
  TVector3 eta = P3.Vect();
  GDouble omega = atan2(y.Dot(eta), P0.Vect().Unit().Dot(eta.Cross(y)));
 

  fillHistogram( PhiT, Phi); 
  fillHistogram( cosT, cosTheta);
  fillHistogram( phiAng, phi);
  fillHistogram( Omega, omega);
  fillHistogram( cosT_m23, (P2+P3).M(), cosTheta);
  fillHistogram( cosT_phi, phi, cosTheta);
  fillHistogram( cosT_Phi, Phi, cosTheta);
}
