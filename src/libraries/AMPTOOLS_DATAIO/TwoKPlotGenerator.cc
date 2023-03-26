#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoKPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

TwoKPlotGenerator::TwoKPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

TwoKPlotGenerator::TwoKPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void TwoKPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( k2PiMass, new Histogram1D( 100, 1.0, 1.04, "M2pi", "Invariant Mass of K^{+} K^{-}") );
  bookHistogram( kPiPCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiPiPlus,  new Histogram1D( 50, -1*PI, PI, "PhiPiPlus",  "#Phi_{K_{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 50, -1*PI, PI, "PhiPiMinus", "#Phi_{K_{-}}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0, 1.00, "t", "-t" ) );

  bookHistogram( kKpLabTheta, new Histogram1D( 100, 0, 20, "KpLabTheta", "Lab #theta of K+ (degrees)") );
  bookHistogram( kKpLabPhi, new Histogram1D( 180, -180, 180, "KpLabPhi", "Lab #phi of K+ (degrees)") );
  bookHistogram( kpKp, new Histogram1D( 300, 0, 9, "pKp", "Lab momentum of K+") );    
  bookHistogram( kKmLabTheta, new Histogram1D( 100, 0, 20, "KmLabTheta", "Lab #theta of K- (degrees)") );
  bookHistogram( kKmLabPhi, new Histogram1D( 180, -180, 180, "KmLabPhi", "Lab #phi of K- (degrees)") );
  bookHistogram( kpKm, new Histogram1D( 300, 0, 9, "pKm", "Lab momentum of K-") );    
  bookHistogram( kPLabTheta, new Histogram1D( 100, 0, 100, "PLabTheta", "Lab #theta of P' (degrees)") );
  bookHistogram( kPLabPhi, new Histogram1D( 180, -180, 180, "PLabPhi", "Lab #phi of P' (degrees)") );  
  bookHistogram( kpP, new Histogram1D( 100, 0, 3, "Pp", "Lab momentum of P'") );      
}

void TwoKPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );

  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  
  // choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
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
  
  fillHistogram( k2PiMass, ( resonance ).M() );
  
  fillHistogram( kPiPCosTheta, cosTheta );

  fillHistogram( kPhiPiPlus,  p1.Phi() );
  fillHistogram( kPhiPiMinus, p2.Phi() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );

  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive


  fillHistogram( kKpLabTheta, 180*p1.Theta()/PI );
  fillHistogram( kKpLabPhi, 180*p1.Phi()/PI );
  fillHistogram( kpKp, p1.P() );
  fillHistogram( kKmLabTheta, 180*p2.Theta()/PI );
  fillHistogram( kKmLabPhi, 180*p2.Phi()/PI );
  fillHistogram( kpKm, p2.P() );
  fillHistogram( kPLabTheta, 180*recoil.Theta()/PI );
  fillHistogram( kPLabPhi, 180*recoil.Phi()/PI );
  fillHistogram( kpP, recoil.P() );


}
