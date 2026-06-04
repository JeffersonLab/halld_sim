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
#include "AMPTOOLS_AMPS/decayAngles.h"

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
   bookHistogram( kLambda, new Histogram1D( 110, 0.0, 1.1, "Lambda", "#lambda_{#omega}" ) );
   bookHistogram( kDalitz, new Histogram2D( 100, -2., 2., 100, -2., 2., "Dalitz", "Dalitz XY" ) );
//   bookHistogram( kThetaDelta, new Histogram1D( 100, 0., PI, "ThetaDelta", "#theta_{#pi^{+}}" ) );
   bookHistogram( kCosThetaDelta, new Histogram1D( 100, -1., 1., "CosThetaDelta", "cos#theta_{#pi^{+}}" ) );
   bookHistogram( kPhiDelta, new Histogram1D( 100, -1*PI, PI, "PhiDelta", "#phi_{#pi^{+}}" ) );
//   bookHistogram( kSinSqThetaDelta, new Histogram1D( 100, 0., 1., "SinSqThetaDelta", "sin^{2}#theta_{#pi^{+}}" ) );
//   bookHistogram( kSin2ThetaDelta, new Histogram1D( 100, -1., 1., "Sin2ThetaDelta", "sin2#theta_{#pi^{+}}" ) );
//   bookHistogram( kCosSqThetaDelta, new Histogram1D( 100, 0., 1., "CosSqThetaDelta", "cos^{2}#theta_{#pi^{+}}" ) );
}

void
OmegaPiPlotGenerator::projectEvent(  Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void
OmegaPiPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector proton = kin->particle( 1 );
   TLorentzVector Xs_pi = kin->particle( 2 );//bachelor pi0
   TLorentzVector omegas_pi = kin->particle( 3 );//omega's pi0
   TLorentzVector rhos_pip = kin->particle( 4 );//pi-
   TLorentzVector rhos_pim = kin->particle( 5 );//pi

   TLorentzVector recoil = proton;
   TLorentzVector proton_pi = proton + Xs_pi;
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
//  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

//  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

   vector< double > upperVertexAngles = getTwoStepAngles( X, omega, rhos_pip, rhos_pim, beam, target, 2, true );
   vector< double > lowerVertexAngles = getOneStepAngles( recoil, proton, beam, target, 2, false );

   GDouble cosTheta = TMath::Cos( upperVertexAngles[0] );
   GDouble Phi = upperVertexAngles[1];
   GDouble cosThetaH = TMath::Cos( upperVertexAngles[2] );
   GDouble PhiH = upperVertexAngles[3];
   GDouble prod_angle = getPhiProd( polAngle, X, beam, target, 2, true );
   GDouble lambda = upperVertexAngles[4];

   GDouble cosThetaDelta = TMath::Cos( lowerVertexAngles[0] );
   GDouble phiDelta = lowerVertexAngles[1];


   // cout << "calls to fillHistogram go here" << endl;
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

   // Dalitz variables
   fillHistogram( kLambda, lambda );
   
   double dalitz_s, dalitz_t, dalitz_u, dalitz_d, dalitz_sc, dalitzx, dalitzy;
   TLorentzVector p2 = omegas_pi;
   TLorentzVector p3 = rhos_pip;
   TLorentzVector p4 = rhos_pim;
   dalitz_s = (p3+p4).M2();//s=M(pip pim)
   dalitz_t = (p2+p3).M2();//s=M(pip pi0)
   dalitz_u = (p2+p4).M2();//s=M(pim pi0)
   dalitz_d = 2*(p2+p3+p4).M()*( (p2+p3+p4).M() - ((2*0.13957018)+0.1349766) );
   dalitz_sc = (1/3.)*( (p2+p3+p4).M2() + ((2*(0.13957018*0.13957018))+(0.1349766*0.1349766)) );
   dalitzx = sqrt(3.)*(dalitz_t - dalitz_u)/dalitz_d;
   dalitzy = 3.*(dalitz_sc - dalitz_s)/dalitz_d;
   fillHistogram( kDalitz, dalitzx, dalitzy );

   // Angles related to lower vertex

   fillHistogram( kCosThetaDelta, cosThetaDelta );
   fillHistogram( kPhiDelta, phiDelta );
}

