#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/EtaPbPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "particleType.h"

/* Constructor to display FitResults */
EtaPbPlotGenerator::EtaPbPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
EtaPbPlotGenerator::EtaPbPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void
EtaPbPlotGenerator::createHistograms( ) {
  // calls to bookHistogram go here 
  bookHistogram( kTheta, new Histogram1D( 50, 0.,5., "Theta", "#theta") );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0.0, 0.2, "t", "-t" ) );
  bookHistogram( kTheta_phi, new Histogram2D( 180, -3.14, 3.14, 100, 0, 5, "Theta_phi", "#theta vs. #phi; #phi; #theta") );
  bookHistogram( kt_phi, new Histogram2D( 180, -3.14, 3.14, 100, 0.0, 0.2, "t_phi", "-t vs. #phi; #phi; -t") );
}

void
EtaPbPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector p1 = kin->particle( 1 );
  TLorentzVector recoil = kin->particle( 2 );

  Double_t kMZ = ParticleMass(Pb208);      //  use mass of Pb as it is in the particle table
  TLorentzVector target  ( 0., 0., 0., kMZ);	
  
  TLorentzVector cm = recoil + p1;
  TLorentzRotation cmBoost( -cm.BoostVector() );
  TLorentzVector p1_cm = cmBoost * p1;
  
  GDouble t = (target - recoil).M2();
  GDouble cosTheta = p1_cm.CosTheta();
  GDouble Theta = 180*acos(cosTheta)/3.14159;
  GDouble phi = p1_cm.Phi();
  if(phi < -1*PI) phi += 2*PI;
  if(phi > PI) phi -= 2*PI;

  cout  << " Plot Generator M=" << beam.M() << " beam="; beam.Print();
  cout << "M=" << p1.M() << " Eta="; p1.Print();
  cout << "M=" << recoil.M() << " -t=" << -t << " recoil="; recoil.Print();

  // calls to fillHistogram go here
  fillHistogram( kTheta, Theta );
  fillHistogram( kPhi, phi );
  fillHistogram( kt, -t );      // fill with -t to make positive
  fillHistogram( kTheta_phi, phi, Theta ); 
  fillHistogram( kt_phi, phi, -t ); 
}
