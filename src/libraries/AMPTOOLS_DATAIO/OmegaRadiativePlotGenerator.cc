#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

#include "AMPTOOLS_DATAIO/OmegaRadiativePlotGenerator.h"
#include "AMPTOOLS_AMPS/wignerD.h"

/* Constructor to display FitResults */
OmegaRadiativePlotGenerator::OmegaRadiativePlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  createHistograms();
}

/* Constructor for event generator (no FitResult) */
OmegaRadiativePlotGenerator::OmegaRadiativePlotGenerator( ) :
PlotGenerator( )
{
  createHistograms();
}

void OmegaRadiativePlotGenerator::createHistograms( ) {

  // calls to bookHistogram go here
  bookHistogram( kOmegaMass, new Histogram1D( 200, 0.72, 0.84, "MOmega", "Invariant Mass;M_{#pi^{0} #gamma}") );
  bookHistogram( kY20, new Histogram1D( 200, 0.72, 0.84, "Y20", "Y^{0}_{2}(#theta); M_{#pi^{0} #gamma}") );
  bookHistogram( kReY21Psi, new Histogram1D( 200, 0.72, 0.84, "ReY21Psi", "Re(Y^{1}_{2}(#theta,#psi)); M_{#pi^{0} #gamma}") );
  bookHistogram( kReY22Psi, new Histogram1D( 200, 0.72, 0.84, "ReY22Psi", "Re(Y^{2}_{2}(#theta,#psi)); M_{#pi^{0} #gamma}") );
  bookHistogram( kReY21Phi, new Histogram1D( 200, 0.72, 0.84, "ReY21Phi", "Re(Y^{1}_{2}(#theta,#phi)); M_{#pi^{0} #gamma}") );
  bookHistogram( kReY22Phi, new Histogram1D( 200, 0.72, 0.84, "ReY22Phi", "Re(Y^{2}_{2}(#theta,#phi)); M_{#pi^{0} #gamma}") );
  bookHistogram( kY40, new Histogram1D( 200, 0.72, 0.84, "Y40", "Y^{0}_{4}(#theta); M_{#pi^{0} #gamma}") );

  bookHistogram( kCosTheta, new Histogram1D( 50, -1, 1, "CosTheta", "cos#theta_{H};cos#theta_{H}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi_{H}; #phi_{H}" ) );
  bookHistogram( kCosThetaPhi, new Histogram2D( 50, -1*PI, PI, 50, -1, 1, "CosThetaPhi", "Helicity Frame; #phi_{H}; cos#theta_{H}" ) );
  bookHistogram( kCosThetaPsi, new Histogram2D( 50, -1*PI, PI, 50, -1, 1, "CosThetaPsi", "cos#theta_{H} vs. #psi; #psi; cos#theta_{H}" ) );
  bookHistogram( kBigPhi, new Histogram1D( 50, -1*PI, PI, "BigPhi", "#Phi; #Phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi; #psi" ) );

  bookHistogram( kt, new Histogram1D( 400, 0, 2.0 , "t", "-t;-t" ) );
  
  bookHistogram( kThetaLabPi0,  new Histogram1D( 40, 0, 20, "ThetaLabPi0",  "#theta lab;#theta_{#pi^{0}}  [deg]" ) );
  bookHistogram( kThetaLabGamma,  new Histogram1D( 40, 0, 20, "ThetaLabGamma",  "#theta lab;#theta_{#gamma}  [deg]" ) );
  bookHistogram( kPThetaLabPi0,  new Histogram2D( 40, 0, 20, 45, 0, 9, "PThetaLabPi0",  ";#theta_{#pi^{0}}  [deg];P_{#pi^{0}}  [GeV]" ) );
  bookHistogram( kPThetaLabGamma,  new Histogram2D( 40, 0, 20, 45, 0, 9, "PThetaLabGamma",  ";#theta_{#gamma}  [deg];P_{#gamma}  [GeV]" ) );
}

void
OmegaRadiativePlotGenerator::projectEvent( Kinematics* pKin, const string& reactionName ){
  
  // obtain the polarzation angle for this event by getting the list of amplitudes
  // associated with this reaction -- we know all are VecRadiative_SDME amplitudes
  // take the tenth argument of the first factor of the first amplitude in the first sum
  // try/catch around stod() to avoid mysterious crash if the conversion from string
  // to double fails
  
  double polAngle = 0;
  string angleString = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).at(10);
  try{ polAngle = stod( angleString ); }
  catch( ... ){
    
    cout << "ERROR:  Unable to get polarization angle from config file because converting\n"
         << "        the amplitude argument '" << angleString << "' to a double failed.\n"
         << "        Check the projectEvent method of the OmegaRadiativePlotGenerator to be sure\n"
         << "        the correct argument is specified.  Also be sure the angle appears\n"
         << "        in the original fit config file as a number and not a parameter.\n"
         << endl;
    assert( false );
  }
    
  // written in the context of omega photoproduction with omega -> gamma  pi0,
  // but these angle definitions in the helicity frame for a two-body decay
  // are fairly standard

  TLorentzVector beam   = pKin->particle( 0 );
  TLorentzVector recoil = pKin->particle( 1 );
  TLorentzVector pi0    = pKin->particle( 2 );
  TLorentzVector gamma  = pKin->particle( 3 );
  
  TLorentzVector omega = pi0 + gamma;
  
  // construct a boost that will take us to the omega rest frame:
  TLorentzRotation omegaRestBoost( -omega.BoostVector() );
  
  TLorentzVector beam_omegaRest = omegaRestBoost * beam;
  TLorentzVector recoil_omegaRest = omegaRestBoost * recoil;
  TLorentzVector gamma_omegaRest = omegaRestBoost * gamma;
  
  // in the helicity frame the z axis is opposite the recoil in the CM
  // frame, but the direction of the recoil in the CM frame is the same
  // as the direction of the recoil in the omega rest frame, so this
  // works fine:
  
  TVector3 zHel = -1*recoil_omegaRest.Vect().Unit();
  
  // the y axis is normal to the production plane production plane, which
  // can be obtained in a variety of ways and is unchanged by boosts
  // to the CM or omega rest frame as all of these boost vectors
  // lie in the production plane
  
  TVector3 yHel = (beam.Vect().Cross(omega.Vect())).Unit();
  
  // and the x axis is constructed so the coordinate system is
  // right handed:
  
  TVector3 xHel = yHel.Cross(zHel);
  
  // a unit vector for the polarization direction in the lab
  TVector3 polUnitVec( cos( polAngle*TMath::DegToRad() ),
                      sin( polAngle*TMath::DegToRad() ), 0.0 );
  
  // now get the angle between the production plane (defined by the
  // y axis in the helicity frame) and the photon polarization vector:
  
  GDouble bigPhi = atan2( yHel.Dot( polUnitVec ),
                          beam.Vect().Unit().Dot( polUnitVec.Cross( yHel ) ) );
  
  // get the angles of the radiated photon in the helicity frame:
  GDouble cosTheta = gamma_omegaRest.Vect().Unit().Dot( zHel );  
  GDouble cosPhi =  yHel.Dot( zHel.Cross( gamma_omegaRest.Vect() ).Unit() );
  GDouble sinPhi = -xHel.Dot( zHel.Cross( gamma_omegaRest.Vect() ).Unit() );
  GDouble phi = atan2( sinPhi, cosPhi );
  
  // this angle is commonly used in old papers for moment calculations
  // it is a difference of angles in two frames, but the form appears
  // in expression for W
  GDouble psi = phi - bigPhi;
  if( psi < -PI ) psi += 2*PI;
  if( psi >  PI ) psi -= 2*PI;
  
  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());
  
  // to use weights in filling, we need a vector of the data that
  // will be filled
  vector< double > omegaMass;
  omegaMass.push_back( omega.M() );
  
  // calls to fillHistogram go here
  fillHistogram( kOmegaMass, omega.M() );
  fillHistogram( kY20, omegaMass, real( Y( 2, 0, cosTheta, psi ) ) );
  fillHistogram( kReY21Psi, omegaMass, real( Y( 2, 1, cosTheta, psi ) ) );
  fillHistogram( kReY22Psi, omegaMass, real( Y( 2, 2, cosTheta, psi ) ) );
  fillHistogram( kReY21Phi, omegaMass, real( Y( 2, 1, cosTheta, phi ) ) );
  fillHistogram( kReY22Phi, omegaMass, real( Y( 2, 2, cosTheta, phi ) ) );
  fillHistogram( kY40, omegaMass, real( Y( 4, 0, cosTheta, psi ) ) );

  fillHistogram( kCosTheta, cosTheta);
  fillHistogram( kPhi, phi );
  fillHistogram( kCosThetaPhi, phi, cosTheta );
  fillHistogram( kCosThetaPsi, psi, cosTheta );
  fillHistogram( kPhi, bigPhi );
  fillHistogram( kPsi, psi );
  
  fillHistogram( kt, -t );      // fill with -t to make positive
  
  fillHistogram( kThetaLabPi0, pi0.Theta()*180./PI );
  fillHistogram( kThetaLabGamma, gamma.Theta()*180./PI );
  fillHistogram( kPThetaLabPi0, pi0.Theta()*180./PI, pi0.P() );
  fillHistogram( kPThetaLabGamma, gamma.Theta()*180./PI, gamma.P() );
}
