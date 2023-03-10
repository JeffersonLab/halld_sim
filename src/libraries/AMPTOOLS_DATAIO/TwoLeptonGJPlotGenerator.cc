#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoLeptonGJPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

TwoLeptonGJPlotGenerator::TwoLeptonGJPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

TwoLeptonGJPlotGenerator::TwoLeptonGJPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void TwoLeptonGJPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( k2LeptonMass, new Histogram1D( 60, 2.9, 3.2, "M2lepton", "Invariant Mass of e^{#minus} e^{#plus}") );
  bookHistogram( kLeptonPCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiLeptonMinus, new Histogram1D( 50, -1*PI, PI, "PhiLeptonMinus", "#Phi_{e_{#minus}}" ) );
  bookHistogram( kPhiLeptonPlus,  new Histogram1D( 50, -1*PI, PI, "PhiLeptonPlus",  "#Phi_{e_{#plus}}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0, 5.00, "t", "-t" ) );
}

void
TwoLeptonGJPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );

  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector beam_res = resonanceBoost * beam;
  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  
  // choose GJ frame: z-axis is defined by the incoming photon in rho rest frame
  TVector3 z = beam_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   (p1_res.Vect()).Dot(x),
                     (p1_res.Vect()).Dot(y),
                     (p1_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  
  GDouble phi = angles.Phi();
  
  TVector3 eps(1.0, 0.0, 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

  // calls to fillHistogram go here
  
  fillHistogram( k2LeptonMass, ( resonance ).M() );
  
  fillHistogram( kLeptonPCosTheta, cosTheta );

  fillHistogram( kPhiLeptonMinus, p1.Phi() );
  fillHistogram( kPhiLeptonPlus,  p2.Phi() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );

  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
