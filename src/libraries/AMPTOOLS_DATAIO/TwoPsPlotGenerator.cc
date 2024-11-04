#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoPsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

TwoPsPlotGenerator::TwoPsPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{

	/*
  vector< string > reactionVec = reactions();
  for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

    // obtain the polarization angle for this reaction by getting the list of amplitudes
    // associated with this reaction -- we know all are SDME amplitudes
    // take the 10th argument of the first factor of the first amplitude in the first sum
    string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(10) };

    // pick out the name of the parameter
    string parName = ampArgument.substr(1,ampArgument.length()-2);

    // here results is of type FitResults which is passed into the constructor of the plot generator
    //cout << "Angle " << results.parValue( parName ) << endl;
    m_reactionAngleMap[*reac] = results.parValue( parName );
  }
	*/

	createHistograms();
}

TwoPsPlotGenerator::TwoPsPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

TwoPsPlotGenerator::~TwoPsPlotGenerator() {
    // No specific cleanup needed; the base class destructor will handle its cleanup
}


void TwoPsPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( k2PsMass, new Histogram1D( 500, 0.5, 2.0, "M2ps", "Invariant Mass of K #pi") );
  bookHistogram( kLambdaKMass, new Histogram1D( 350, 1, 4.5, "MLambdaK", "Invariant Mass of #Lambda K") );
  bookHistogram( kLambdaPiMass, new Histogram1D( 350, 1, 4.5, "MLambdaPi", "Invariant Mass of #Lambda #pi") );
  bookHistogram( kPiCosTheta, new Histogram1D( 200, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiK,  new Histogram1D( 180, -1*PI, PI, "PhiK",  "#Phi_{K}" ) );
  bookHistogram( kPhiPi, new Histogram1D( 180, -1*PI, PI, "PhiPi", "#Phi_{#pi}" ) );
  bookHistogram( kPhiLambda,  new Histogram1D( 180, -1*PI, PI, "PhiLambda", "#Phi_{#Lambda}" ) );
  
  bookHistogram( kThetaK,  new Histogram1D( 200, 0, 20, "ThetaK",  "#Theta_{K}" ) );
  bookHistogram( kThetaPi, new Histogram1D( 200, 0, 20, "ThetaPi", "#Theta_{#pi}" ) );
  bookHistogram( kThetaLambda, new Histogram1D( 200, 50, 90, "ThetaLambda", "#Theta_{#Lambda}" ) );
  
  bookHistogram( kMomK,  new Histogram1D( 180, 0, 9, "MomK",  "p_{K}" ) );
  bookHistogram( kMomPi, new Histogram1D( 180, 0, 9, "MomPi", "p_{#pi}" ) );
  bookHistogram( kMomLambda, new Histogram1D( 180, 0, 3, "MomLambda", "p_{#Lambda}" ) );
  
  bookHistogram( kPhi, new Histogram1D( 180, -1*PI, PI, "Phi", "#Phi Polarization (Lab Frame)" ) );
  bookHistogram( kphi, new Histogram1D( 180, -1*PI, PI, "phi", "#phi (Helicity Frame)" ) );
  bookHistogram( kPsi, new Histogram1D( 180, -1*PI, PI, "psi", "#psi (#phi-#Phi)" ) );
  bookHistogram( kt, new Histogram1D( 500, 0, 1.00, "t", "-t" ) );
}

void TwoPsPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used

  projectEvent( kin, "" );
}

void TwoPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  double polAngle = 0.0; //m_reactionAngleMap[ reactionName ];
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 ); // Lambda
  TLorentzVector p1 = kin->particle( 2 ); // K
  TLorentzVector p2 = kin->particle( 3 ); // Pi

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
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

  // calls to fillHistogram go here
  
  fillHistogram( k2PsMass, ( resonance ).M() );
  fillHistogram( kLambdaKMass, ( recoil+p1 ).M() );
  fillHistogram( kLambdaPiMass, ( recoil+p2 ).M() );
  fillHistogram( kPiCosTheta, cosTheta );
  fillHistogram( kPhiK,  p1.Phi() );
  fillHistogram( kPhiPi, p2.Phi() );
  fillHistogram( kPhiLambda, recoil.Phi() );
  fillHistogram( kThetaK,  p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPi, p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaLambda, recoil.Theta()*TMath::RadToDeg() );
  fillHistogram( kMomK,  p1.P() );
  fillHistogram( kMomPi, p2.P() );
  fillHistogram( kMomLambda, recoil.P() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );
  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
