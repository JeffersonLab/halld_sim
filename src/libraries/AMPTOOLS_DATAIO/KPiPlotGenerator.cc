#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/KPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

KPiPlotGenerator::KPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  vector< string > reactionVec = reactions();
  for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

    // obtain the polarization angle for this reaction by getting the list of amplitudes
    // associated with this reaction -- we know all are SDME amplitudes
    // take the 10th argument of the first factor of the first amplitude in the first sum
    string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(10) };

    // pick out the name of the parameter
    string parName = ampArgument.substr(1, ampArgument.length()-2);

    m_reactionAngleMap[*reac] = results.parValue( parName );
  }

  createHistograms();
}

KPiPlotGenerator::KPiPlotGenerator( ) :
PlotGenerator( )
{
  createHistograms();
}

void KPiPlotGenerator::createHistograms() {

  // reaction:  Beam(0)  Lambda(1)  KShort(2)  Pi+(3)
  // resonance: K*(892) --> KShort + Pi+

  bookHistogram( kKPiMass,      new Histogram1D( 200, 0.6,  2.5,  "MKPi",      "Invariant Mass of K^{0}_{S}#pi^{+}" ) );
  bookHistogram( kLambKMass,    new Histogram1D( 200, 1.4,  4.0,  "MLambK",    "Invariant Mass of #Lambda K^{0}_{S}" ) );
  bookHistogram( kLambPiMass,   new Histogram1D( 200, 1.1,  4.0,  "MLambPi",   "Invariant Mass of #Lambda#pi^{+}" ) );
  bookHistogram( kKCosTheta,    new Histogram1D(  36, -1.,  1.,   "cosTheta",  "cos( #theta ) of K_{s}^{+}" ) );
  bookHistogram( kThetaK,       new Histogram1D( 200,  0,   20,   "ThetaK",    "#Theta_{K^{0}_{S}}" ) );
  bookHistogram( kThetaPi,      new Histogram1D( 200,  0,   20,   "ThetaPi",   "#Theta_{#pi^{+}}" ) );
  bookHistogram( kThetaLamb,    new Histogram1D( 200, 50,   90,   "ThetaLamb", "#Theta_{#Lambda}" ) );
  bookHistogram( kMomK,         new Histogram1D( 180,  0,    9,   "MomK",      "p_{K^{0}_{S}}" ) );
  bookHistogram( kMomPi,        new Histogram1D( 180,  0,    9,   "MomPi",     "p_{#pi^{+}}" ) );
  bookHistogram( kMomLamb,      new Histogram1D( 180,  0,    3,   "MomLamb",   "p_{#Lambda}" ) );
  bookHistogram( kPhiK,         new Histogram1D( 180, -PI,  PI,   "PhiK",      "#phi_{K^{0}_{S}}" ) );
  bookHistogram( kPhiPi,        new Histogram1D( 180, -PI,  PI,   "PhiPi",     "#phi_{#pi^{+}}" ) );
  bookHistogram( kPhiLamb,      new Histogram1D( 180, -PI,  PI,   "PhiLamb",   "#phi_{#Lambda}" ) );
  bookHistogram( kPhi,          new Histogram1D(  36, -PI,  PI,   "Phi",       "#Phi" ) );
  bookHistogram( kphi,          new Histogram1D(  36, -PI,  PI,   "phi",       "#phi" ) );
  bookHistogram( kPsi,          new Histogram1D(  36, -PI,  PI,   "psi",       "#psi" ) );
  bookHistogram( kt,            new Histogram1D( 500,  0,   1.0,  "t",         "-t" ) );
}

void KPiPlotGenerator::projectEvent( Kinematics* kin ){

  // backwards-compatible with older AmpTools versions
  projectEvent( kin, "" );
}

void KPiPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  double polAngle = m_reactionAngleMap[ reactionName ];

  // particle ordering matches config:  Beam(0)  Lambda(1)  KShort(2)  Pi+(3)
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );  // Lambda
  TLorentzVector p1     = kin->particle( 2 );  // KShort
  TLorentzVector p2     = kin->particle( 3 );  // Pi+

  TLorentzVector resonance = p1 + p2;          // K*(892) candidate
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res     = resonanceBoost * p1;

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

  // helicity frame: z-axis opposite recoil Lambda in K* rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();

  TVector3 angles( (p1_res.Vect()).Dot(x),
                   (p1_res.Vect()).Dot(y),
                   (p1_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  GDouble phi      = angles.Phi();

  TVector3 eps( cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0 );
  GDouble Phi = atan2( y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)) );

  GDouble psi = phi - Phi;
  if( psi < -PI ) psi += 2*PI;
  if( psi >  PI ) psi -= 2*PI;

  // invariant t
  GDouble t = -2.0 * recoil.M() * ( recoil.E() - recoil.M() );

  // fill histograms
  fillHistogram( kKPiMass,       resonance.M() );
  fillHistogram( kLambKMass,     ( recoil + p1 ).M() );   // Lambda + KShort
  fillHistogram( kLambPiMass,    ( recoil + p2 ).M() );   // Lambda + Pi+
  fillHistogram( kKCosTheta,     cosTheta );
  fillHistogram( kThetaK,        p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPi,       p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaLamb,     recoil.Theta()*TMath::RadToDeg() );
  fillHistogram( kMomK,          p1.P() );
  fillHistogram( kMomPi,         p2.P() );
  fillHistogram( kMomLamb,       recoil.P() );
  fillHistogram( kPhiK,          p1.Phi() );
  fillHistogram( kPhiPi,         p2.Phi() );
  fillHistogram( kPhiLamb,       recoil.Phi() );
  fillHistogram( kPhi,           Phi );
  fillHistogram( kphi,           phi );
  fillHistogram( kPsi,           psi );
  fillHistogram( kt,             -t );    // fill with -t to make positive
}
