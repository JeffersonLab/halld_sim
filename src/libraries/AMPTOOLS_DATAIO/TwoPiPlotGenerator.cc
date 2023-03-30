#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

TwoPiPlotGenerator::TwoPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{

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

	createHistograms();
}

TwoPiPlotGenerator::TwoPiPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void TwoPiPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( k2PiMass, new Histogram1D( 500, 0.5, 1.0, "M2pi", "Invariant Mass of #pi^{+} #pi^{-}") );
  bookHistogram( kPPipMass, new Histogram1D( 350, 1, 4.5, "Mppip", "Invariant Mass of p#pi^{+}") );
  bookHistogram( kPPimMass, new Histogram1D( 350, 1, 4.5, "Mppim", "Invariant Mass of p#pi^{-}") );
  bookHistogram( kPiPCosTheta, new Histogram1D( 200, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiPiPlus,  new Histogram1D( 180, -1*PI, PI, "PhiPiPlus",  "#Phi_{#pi_{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 180, -1*PI, PI, "PhiPiMinus", "#Phi_{#pi_{-}}" ) );
  bookHistogram( kPhiProton,  new Histogram1D( 180, -1*PI, PI, "PhiProton", "#Phi_{p}" ) );
  bookHistogram( kThetaPiPlus,  new Histogram1D( 200, 0, 20, "ThetaPiPlus",  "#Theta_{#pi_{+}}" ) );
  bookHistogram( kThetaPiMinus, new Histogram1D( 200, 0, 20, "ThetaPiMinus", "#Theta_{#pi_{-}}" ) );
  bookHistogram( kThetaProton, new Histogram1D( 200, 50, 90, "ThetaProton", "#Theta_{p}" ) );
  bookHistogram( kMomPiPlus,  new Histogram1D( 180, 0, 9, "MomPiPlus",  "p_{#pi_{+}}" ) );
  bookHistogram( kMomPiMinus, new Histogram1D( 180, 0, 9, "MomPiMinus", "p_{#pi_{-}}" ) );
  bookHistogram( kMomProton, new Histogram1D( 180, 0, 3, "MomProton", "p_{p}" ) );
  bookHistogram( kPhi, new Histogram1D( 180, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 180, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 180, -1*PI, PI, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 500, 0, 1.00, "t", "-t" ) );
}

void TwoPiPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used

  projectEvent( kin, "" );
}

void TwoPiPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  double polAngle = m_reactionAngleMap[ reactionName ];
  
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
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

  // calls to fillHistogram go here
  
  fillHistogram( k2PiMass, ( resonance ).M() );
  fillHistogram( kPPipMass, ( recoil+p1 ).M() );
  fillHistogram( kPPimMass, ( recoil+p2 ).M() );
  
  fillHistogram( kPiPCosTheta, cosTheta );

  fillHistogram( kPhiPiPlus,  p1.Phi() );
  fillHistogram( kPhiPiMinus, p2.Phi() );
  fillHistogram( kPhiProton, recoil.Phi() );
  fillHistogram( kThetaPiPlus,  p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPiMinus, p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaProton, recoil.Theta()*TMath::RadToDeg() );
  fillHistogram( kMomPiPlus,  p1.P() );
  fillHistogram( kMomPiMinus, p2.P() );
  fillHistogram( kMomProton, recoil.P() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );

  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
