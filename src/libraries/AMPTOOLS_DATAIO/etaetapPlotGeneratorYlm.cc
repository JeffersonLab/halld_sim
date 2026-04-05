#include "AMPTOOLS_DATAIO/etaetapPlotGeneratorYlm.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Ylm.h"
#include "AMPTOOLS_AMPS/Ylm_reim.h"


etaetapPlotGeneratorYlm::etaetapPlotGeneratorYlm( const FitResults& results ) :
PlotGenerator( results )
{

  bookHistogram( khm12, new Histogram1D( 200, 1.0, 5.0, "hm12", "Mass( 1 2 )" ) );
  bookHistogram( khm13, new Histogram1D( 200, 1.0, 5.0, "hm13", "Mass( 1 3 )" ) );
  bookHistogram( khm23, new Histogram1D( 20, 1.4, 3.4, "hm23", "M(#eta#eta') [GeV/c^{2}]" ) );
  bookHistogram( khm1, new Histogram1D( 100, 0.5, 1.5, "hm1", "Mass( 1 )" ) );
  bookHistogram( khm2, new Histogram1D( 100, 0.0, 1.0, "hm2", "Mass( 2 )" ) );
  bookHistogram( khm3, new Histogram1D( 100, 0.0, 1.0, "hm3", "Mass( 3 )" ) );
  bookHistogram( kdltz, new Histogram2D( 80, 0.0, 25.0, 80, 0.0, 9.0, "dltz", "Dalitz Plot" ) );
  bookHistogram( cosT, new Histogram1D( 55, -1.1, 1.1, "cosT", "CosTheta") );
  bookHistogram( phiAng, new Histogram1D(50, -3.2, 3.2, "phiAng", "#phi") );
  bookHistogram( cosT_m23, new Histogram2D(100, 0.7, 2.7, 100, -1.0, 1.0, "cosT_m23", "cos(#theta) vs. Mass(#eta #eta')" ) );
  bookHistogram( cosT_phi, new Histogram2D(40, -3.2, 3.2, 100, -1.0, 1.0, "cosT_phi", "cos(#theta) vs. #phi" ) );
}

void etaetapPlotGeneratorYlm::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void etaetapPlotGeneratorYlm::projectEvent( Kinematics* kin, const string& reactionName ){

  TLorentzVector P0 = kin->particle(0);// beam photon
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
  
  TVector3 eta = P3.Vect();
 

  fillHistogram( cosT, cosTheta);
  fillHistogram( phiAng, phi);
  fillHistogram( cosT_m23, (P2+P3).M(), cosTheta);
  fillHistogram( cosT_phi, phi, cosTheta);

}
