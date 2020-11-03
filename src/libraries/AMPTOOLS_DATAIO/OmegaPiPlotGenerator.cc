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
OmegaPiPlotGenerator::OmegaPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
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
  
   bookHistogram( kOmegaPiMass, new Histogram1D( 240, 0.8, 3., "MOmegaPi", "Invariant Mass of #omega #pi^{0}") );
   bookHistogram( kOmegaMass, new Histogram1D( 120, 0.4, 1., "MOmega", "Invariant Mass of #pi^{+} #pi^{-} #pi^{0}") );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta;cos#theta" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi; #Phi[rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta_H", "cos#theta_H;cos#theta_H" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -1*PI, PI, "Phi_H", "#Phi_H; #Phi_H[rad.]" ) );
   bookHistogram( kProd_Ang, new Histogram1D( 50, -1*PI, PI, "Prod_Ang", "Prod_Ang; Prod_Ang[rad.]" ) );
   bookHistogram( kt, new Histogram1D( 100, 0, 2.0 , "t", "-t" ) );
   
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

   for(uint i=6; i<kin->particleList().size(); i++) recoil += kin->particle(i);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TLorentzVector rho = rhos_pip + rhos_pim;
  TLorentzVector omega = rho + omegas_pi;
  TLorentzVector X = omega + Xs_pi;

  TLorentzVector target(0,0,0,0.938);
  double polAngle = 0.0;
  //cout << "Beam Polarization Angle Set to " << polAngle << endl;

  double b1_mass = X.M();
  double omega_mass = omega.M();
  double Mandt = fabs((target-recoil).M2());
  
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
   fillHistogram( kOmegaMass, omega_mass );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, Phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, PhiH );
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kt, Mandt );

}
