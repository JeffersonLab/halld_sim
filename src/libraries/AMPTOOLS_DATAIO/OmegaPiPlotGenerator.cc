//#include <ctime>
//#include <stdlib.h>
//#include <stdio.h>

//#include <cassert>
//#include <iostream>
//#include <string>
//#include <sstream>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/omegapiAngles.h"

//#include <cmath>
//#include <complex>
//#include <vector>
//#include "TMath.h"

#include "AMPTOOLS_DATAIO/OmegaPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

/* Constructor to display FitResults */
OmegaPiPlotGenerator::OmegaPiPlotGenerator( const FitResults& results, Option opt ) :
PlotGenerator( results, opt )
{
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
OmegaPiPlotGenerator::OmegaPiPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void OmegaPiPlotGenerator::createHistograms( ) {
  cout << " calls to bookHistogram go here" << endl;
  
   bookHistogram( kOmegaPiMass, new Histogram1D( 200, 0.6, 2., "MOmegaPi", "Invariant Mass of #omega #pi") );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#phi[rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta_H", "cos#theta_H" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -1*PI, PI, "Phi_H", "#phi_H[rad.]" ) );
   bookHistogram( kProd_Ang, new Histogram1D( 50, -1*PI, PI, "Prod_Ang", "Prod_Ang[rad.]" ) );
   bookHistogram( kt, new Histogram1D( 100, 0, 2.0 , "t", "-t" ) );
   bookHistogram( kRecoilMass, new Histogram1D( 100, 0.9, 1.9 , "MRecoil", "Invariant Mass of Recoil" ) );
   bookHistogram( kTwoPiMass, new Histogram1D( 100, 0.25, 1.75, "MTwoPi", "Invariant Mass of 2 #pi" ) ); 
   bookHistogram( kProtonPiMass, new Histogram1D( 100, 0.9, 2.9, "MProtonPi", "Invariant Mass of proton and bachelor pion" ) );
   bookHistogram( kRecoilPiMass, new Histogram1D( 100, 0.9, 2.9, "MRecoilPi", "Invariant Mass of recoil and bachelor pion" ) );
  
}

void
OmegaPiPlotGenerator::projectEvent( Kinematics* kin ){

   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector Xs_pi = kin->particle( 2 );//bachelor pi0
   TLorentzVector omegas_pi = kin->particle( 3 );//omega's pi0
   TLorentzVector rhos_pip = kin->particle( 4 );//pi-
   TLorentzVector rhos_pim = kin->particle( 5 );//pi+

   TLorentzVector proton_pi = recoil + Xs_pi;
   TLorentzVector recoil_pi = proton_pi;
   TLorentzVector two_pi = kin->particle( 2 );
   for(uint i=6; i<kin->particleList().size(); i++) {
	recoil += kin->particle(i);
	recoil_pi += kin->particle(i);
	two_pi += kin->particle(i);
   }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TLorentzVector rho = rhos_pip + rhos_pim;
  TLorentzVector omega = rho + omegas_pi;
  TLorentzVector X = omega + Xs_pi;

  TLorentzVector target(0,0,0,0.938);
  double polAngle = 0.0;
  //cout << "Beam Polarization Angle Set to " << polAngle << endl;

  double b1_mass = X.M();
  double Mandt = fabs((target-recoil).M2());
  double recoil_mass = recoil.M();  

    //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;
 
  //Calculate decay angles in helicity frame
  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

   GDouble cosTheta = TMath::Cos(locthetaphi[0]);
   GDouble Phi = locthetaphi[1];
   GDouble cosThetaH = TMath::Cos(locthetaphih[0]);
   GDouble PhiH = locthetaphih[1];
   GDouble prod_angle = locthetaphi[2];

   //cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kOmegaPiMass, b1_mass );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, Phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, PhiH );
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kt, Mandt );
   fillHistogram( kRecoilMass, recoil_mass );
   fillHistogram( kTwoPiMass, two_pi.M() );
   fillHistogram( kProtonPiMass, proton_pi.M() );
   fillHistogram( kRecoilPiMass, recoil_pi.M() );

}

