#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoPiDeltaPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

TwoPiDeltaPlotGenerator::TwoPiDeltaPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{

  vector< string > reactionVec = reactions();
  for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

    // obtain the polarization angle for this reaction by getting the list of amplitudes
    // associated with this reaction -- we know all are SDME amplitudes
    // // take the 1th argument of the first factor of the first amplitude in the first sum
    // string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(1) };

    // pick out the name of the parameter
    // string parName = ampArgument.substr(1,ampArgument.length()-2);

    // // here results is of type FitResults which is passed into the constructor of the plot generator
    // //cout << "Angle " << results.parValue( parName ) << endl;
    // m_reactionAngleMap[*reac] = results.parValue( parName );
  
    double polAngle = 0;
    string angleString = cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(2);
   
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
  m_reactionAngleMap[*reac] = polAngle;
  

  }



for(const auto& entry : m_reactionAngleMap){
  cout << "DEBUG: Reaction in map: " << entry.first 
       << "  with angle: " << entry.second << endl;
}

	createHistograms();
}

TwoPiDeltaPlotGenerator::TwoPiDeltaPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void TwoPiDeltaPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  //signal mass
  bookHistogram( k2PiMass, new Histogram1D( 50, 0.6, 0.95, "M2pi", "M(#pi^{0} #pi^{-}) (GeV)") );
  bookHistogram( kPPipMass, new Histogram1D( 50, 1.1, 1.4, "Mppip", "M(p#pi^{+} )(GeV)") );
   //different combinations of invariant masses
  bookHistogram( kPipPimMass, new Histogram1D( 100, 0.2, 2.5, "Mpippim", "M(#pi^{+}#pi^{-} )(GeV)") );
  bookHistogram( kPPi0Mass, new Histogram1D( 100, 1, 4.5, "Mppi0", "M(p#pi^{0}) (GeV)") );
  bookHistogram( kPipPi0Mass, new Histogram1D( 100, 0.1, 2.5, "Mpippi0", "M(#pi^{+}#pi^{0}) (GeV)") );
  
  bookHistogram( kPPimMass, new Histogram1D( 100, 1, 4.5, "Mppim", "M( p#pi^{-}) (GeV)") );
  bookHistogram( kPi0PimMass, new Histogram1D( 100, 0.1, 1.5, "Mpi0pim", "M(#pi^{0}#pi^{-}) (GeV)") );
  bookHistogram( kPPipPi0Mass, new Histogram1D( 100, 1, 4.5, "Mppippi0", "M( p#pi^{+}#pi^{0}) (GeV)") );
  bookHistogram( kPPimPi0Mass, new Histogram1D( 100, 2.5, 4.5, "Mppimpi0", "M( p#pi^{-}#pi^{0}) (GeV)") );
  bookHistogram( kPipPi0PimMass, new Histogram1D( 100, 0.5, 4.5, "Mpippi0pim", "M( #pi^{-}#pi^{0}#pi^{+}) (GeV)") ); 
  bookHistogram( kPPipPimMass, new Histogram1D( 100, 1, 4.5, "Mppippim", "M( p#pi^{-}#pi^{+}) (GeV)") );


 
  bookHistogram( kPhiPiPlus,  new Histogram1D( 180, -180, 180, "PhiPiPlus",  "#Phi_{#pi^{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 180, -180, 180, "PhiPiMinus", "#Phi_{#pi^{-}}" ) );
  bookHistogram( kPhiProton,  new Histogram1D( 180, -180, 180, "PhiProton", "#Phi_{p}" ) );

  bookHistogram( kThetaPiPlus,  new Histogram1D( 200, 0, 20, "ThetaPiPlus",  "#Theta_{#pi^{+}}" ) );
  bookHistogram( kThetaPiMinus, new Histogram1D( 200, 0, 20, "ThetaPiMinus", "#Theta_{#pi^{-}}" ) );
  bookHistogram( kThetaDelta, new Histogram1D( 200, 50, 90, "ThetaDelta", "#Theta_{#Delta}" ) );

  
  

  // longitudinal momenta 
  bookHistogram( klongMomPiPlus,  new Histogram1D( 60, 0, 0.5, "longMomPiPlus",  "p_{z,#pi^{+}}" ) );
  bookHistogram( klongMomPiMinus, new Histogram1D( 60, 0, 9, "longMomPiMinus", "p_{z,#pi^{-}}" ) );
  bookHistogram( klongMomPi0, new Histogram1D( 60, 0, 9, "longMomPi0", "p_{z,#pi^{0}}" ) );
  bookHistogram( klongMomProton, new Histogram1D( 60, 0, 1.5, "longMomProton", "p_{z,p}" ) );

  // momenta 
  bookHistogram( kMomPiPlus,  new Histogram1D( 60, 0, 0.7, "MomPiPlus",  "p_{#pi^{+}}" ) );
  bookHistogram( kMomPiMinus, new Histogram1D( 60, 0, 9, "MomPiMinus", "p_{#pi^{-}}" ) );
  bookHistogram( kMomPi0, new Histogram1D( 60, 0, 9, "MomPi0", "p_{#pi^{0}}" ) );
  bookHistogram( kMomProton, new Histogram1D( 180, 0, 1, "MomProton", "p_{p}" ) );
  
  // Angles in different frames Phi
  bookHistogram( kphi_PiMinus_hel_rho, new Histogram1D( 60, -180, 180, "phi_PiMinus_hel_rho", " #phi_{#pi^{-}} in hel (restframe #rho) " ) );
  bookHistogram( kphi_PiMinus_GJ_rho, new Histogram1D( 60, -180, 180, "phi_PiMinus_GJ_rho", "#phi_{#pi^{-}} in GJ (restframe #rho)" ) );

  bookHistogram( kphi_PiPlus_hel_Delta, new Histogram1D( 60, -180, 180, "phi_PiPlus_hel_Delta", " #phi_{#pi^{+}} in hel (restframe #Delta)" ) );
  bookHistogram( kphi_PiPlus_GJ_Delta, new Histogram1D( 60, -180, 180, "phi_PiPlus_GJ_Delta", "#phi_{#pi^{+}} in GJ (restframe #Delta)" ) );

  bookHistogram( kphi_Proton_hel_Delta, new Histogram1D( 60, -180, 180, "phi_Proton_hel_Delta", " #phi_{p} in hel (restframe #Delta)" ) );
  bookHistogram( kphi_Proton_GJ_Delta, new Histogram1D( 60, -180, 180, "phi_Proton_GJ_Delta", "#phi_{p} in GJ (restframe #Delta)" ) );

  bookHistogram( kCosTheta_PiMinus_hel_rho, new Histogram1D( 60, -1., 1., "cosTheta_PiMinus_hel_rho", "cos( #theta_{#pi^{-}} ) in hel (restframe #rho)") );
  bookHistogram( kCosTheta_PiMinus_GJ_rho, new Histogram1D( 60, -1., 1., "cosTheta_PiMinus_GJ_rho", "cos( #theta_{#pi^{-}} ) in GJ (restframe #rho)") );

  bookHistogram( kCosTheta_PiPlus_hel_Delta, new Histogram1D( 60, -1., 1., "cosTheta_PiPlus_hel_Delta", "cos( #theta_{#pi^{+}} ) in hel (restframe #Delta)") );
  bookHistogram( kCosTheta_PiPlus_GJ_Delta, new Histogram1D( 60, -1., 1., "cosTheta_PiPlus_GJ_Delta", "cos( #theta_{#pi^{+}} ) in GJ (restframe #Delta)") );

  bookHistogram( kCosTheta_Proton_hel_Delta, new Histogram1D( 60, -1., 1., "cosTheta_Proton_hel_Delta", "cos( #theta_{p} ) in hel (restframe #Delta)") ); 
  bookHistogram( kCosTheta_Proton_GJ_Delta, new Histogram1D( 60, -1., 1., "cosTheta_Proton_GJ_Delta", "cos( #theta_{p} ) in GJ (restframe #Delta)") );
  

  bookHistogram( kPhi, new Histogram1D( 50, -180, 180, "Phi", "#Phi" ) );
  bookHistogram( kPsi, new Histogram1D( 20, -180, 180, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 1000, 0, 1.4, "-t", "-t (GeV^2)" ) );
  
}

void TwoPiDeltaPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used

  projectEvent( kin, "" );
}

void TwoPiDeltaPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  double polAngle = m_reactionAngleMap[ reactionName ];
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector proton = kin->particle( 1 ); //proton
  TLorentzVector p1 = kin->particle( 2 ); //pip 
  TLorentzVector p2 = kin->particle( 3 ); //pim
  TLorentzVector p3 = kin->particle( 4 ); //pi0 
  TLorentzVector target ( 0, 0, 0, 0.9382720813);

  TLorentzVector recoil = proton + p1;
  TLorentzVector rho = p2 + p3;
  TLorentzVector Delta = proton + p1;
  TLorentzVector CM_motion_lab = beam + target;
  TLorentzVector piminus = p2;

  TLorentzRotation rhoBoost( -rho.BoostVector() );
  TLorentzRotation CMBoost( -CM_motion_lab.BoostVector() );
  TLorentzRotation DeltaBoost( -Delta.BoostVector() );


// choose helicity frame: z-axis opposite recoil Delta in rho rest frame
//bost in rho restframe
  TLorentzVector beam_resrho = rhoBoost * beam;
  TLorentzVector Delta_resrho = rhoBoost * Delta;
  TLorentzVector piminus_resrho = rhoBoost * piminus;
  TLorentzVector proton_resrho = rhoBoost * proton;
//boost in CM frame
  TLorentzVector beam_cm = CMBoost * beam;
  TLorentzVector Delta_cm = CMBoost * Delta;
  TLorentzVector piminus_cm = CMBoost * piminus;
  TLorentzVector rho_cm = CMBoost * rho;
  
//boost in recoil frame
  TLorentzVector target_resDelta = DeltaBoost * target;
  TLorentzVector proton_resDelta = DeltaBoost * proton;
  TLorentzVector piplus_resDelta = DeltaBoost * p1;
  TLorentzVector rho_resDelta = DeltaBoost * rho;
  
  // Define coordinate systems: 
  // y equal in all frames
  TVector3 y = (beam.Vect().Unit().Cross(-Delta.Vect().Unit())).Unit();

  //rho system
  
  //Helicity frame
  TVector3 z_hel_rho = -1. * Delta_resrho.Vect().Unit();
  TVector3 x_hel_rho = y.Cross(z_hel_rho).Unit();

  // GJ frame

  TVector3 z_GJ_rho = beam_resrho.Vect().Unit();
  TVector3 x_GJ_rho = y.Cross(z_GJ_rho).Unit();

  //Delta system

  // Helicity frame
  TVector3 z_hel_Delta = -rho_resDelta.Vect().Unit();
	TVector3 x_hel_Delta = (y.Cross(z_hel_Delta)).Unit();

  // GJ frame
  TVector3 z_GJ_Delta = target_resDelta.Vect().Unit();
  TVector3 x_GJ_Delta = y.Cross(z_GJ_Delta).Unit();

  //Angle calculation
  TVector3 angles_piminus_hel_rho( (piminus_resrho.Vect()).Dot(x_hel_rho),
		   (piminus_resrho.Vect()).Dot(y),
		   (piminus_resrho.Vect()).Dot(z_hel_rho) );
  
  TVector3 angles_piminus_GJ_rho( (piminus_resrho.Vect()).Dot(x_GJ_rho),
		   (piminus_resrho.Vect()).Dot(y),
		   (piminus_resrho.Vect()).Dot(z_GJ_rho) );


  TVector3 angles_proton_hel_Delta( (proton_resDelta.Vect()).Dot(x_hel_Delta),
		   (proton_resDelta.Vect()).Dot(y),
		   (proton_resDelta.Vect()).Dot(z_hel_Delta) );
  TVector3 angles_piplus_hel_Delta( (piplus_resDelta.Vect()).Dot(x_hel_Delta),
		   (piplus_resDelta.Vect()).Dot(y),
		   (piplus_resDelta.Vect()).Dot(z_hel_Delta) );
  
  TVector3 angles_proton_GJ_Delta( (proton_resDelta.Vect()).Dot(x_GJ_Delta),
		   (proton_resDelta.Vect()).Dot(y),
		   (proton_resDelta.Vect()).Dot(z_GJ_Delta) );

  TVector3 angles_piplus_GJ_Delta( (piplus_resDelta.Vect()).Dot(x_GJ_Delta),
  (piplus_resDelta.Vect()).Dot(y),
  (piplus_resDelta.Vect()).Dot(z_GJ_Delta) );
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = angles_piminus_GJ_rho.Phi() - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = -(beam-p3-p2).M2();

  // calls to fillHistogram go here
  fillHistogram( k2PiMass, ( rho ).M() );
  fillHistogram( kPPipMass, ( proton+p1 ).M() );
  fillHistogram( kPPimMass, ( proton+p2 ).M() );
  fillHistogram( kPi0PimMass, ( p3+p2 ).M() );
  
 
  fillHistogram( kPipPimMass,   (p1+p2 ).M());
  fillHistogram( kPPi0Mass,   (proton+p3 ).M());
  fillHistogram( kPipPi0Mass,   (p1+p3 ).M());

  ;
  fillHistogram( klongMomPiPlus,  p1.Pz() );
  fillHistogram( klongMomPiMinus, p2.Pz() );
  fillHistogram( klongMomPi0, p3.Pz() );
  fillHistogram( klongMomProton, proton.Pz() );

  fillHistogram( kMomPiPlus,  p1.P() );
  fillHistogram( kMomPiMinus, p2.P() );
  fillHistogram( kMomPi0, p3.P() );
  fillHistogram( kMomProton, proton.P() );
  fillHistogram( kPPipPi0Mass, ( proton + p1 + p3 ).M() );
  fillHistogram( kPipPi0PimMass, ( p1 + p3 + p2 ).M() );
  fillHistogram( kPPipPimMass, ( proton + p1 + p2 ).M() );
  fillHistogram( kPPimPi0Mass, ( proton + p2 + p3 ).M() );
  

  //angles for fits

  // phi angles in different frames
  fillHistogram( kphi_PiMinus_hel_rho, angles_piminus_hel_rho.Phi()*TMath::RadToDeg() );
  fillHistogram( kphi_PiMinus_GJ_rho, angles_piminus_GJ_rho.Phi()*TMath::RadToDeg());

  fillHistogram( kphi_PiPlus_hel_Delta, angles_piplus_hel_Delta.Phi()*TMath::RadToDeg());
  fillHistogram( kphi_PiPlus_GJ_Delta, angles_piplus_GJ_Delta.Phi()*TMath::RadToDeg());

  fillHistogram( kphi_Proton_hel_Delta, angles_proton_hel_Delta.Phi()*TMath::RadToDeg());
  fillHistogram( kphi_Proton_GJ_Delta, angles_proton_GJ_Delta.Phi()*TMath::RadToDeg()); 
  
  // cosTheta angles in different frames
  fillHistogram( kCosTheta_PiMinus_hel_rho, angles_piminus_hel_rho.CosTheta() );
  fillHistogram( kCosTheta_PiMinus_GJ_rho, angles_piminus_GJ_rho.CosTheta() );

  fillHistogram( kCosTheta_PiPlus_hel_Delta, angles_piplus_hel_Delta.CosTheta() );
  fillHistogram( kCosTheta_PiPlus_GJ_Delta, angles_piplus_GJ_Delta.CosTheta() );
  
  fillHistogram( kCosTheta_Proton_hel_Delta, angles_proton_hel_Delta.CosTheta() ); 
  fillHistogram( kCosTheta_Proton_GJ_Delta, angles_proton_GJ_Delta.CosTheta() );

  //Different angles in lab frame
  fillHistogram( kThetaPiPlus,  p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPiMinus, p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaDelta, recoil.Theta()*TMath::RadToDeg() );
  fillHistogram( kPhiPiPlus,  p1.Phi()*TMath::RadToDeg() );
  fillHistogram( kPhiPiMinus, p2.Phi()*TMath::RadToDeg() );
  fillHistogram( kPhiProton, proton.Phi()*TMath::RadToDeg() );

  //Polarisation angles
  fillHistogram( kPhi, Phi*TMath::RadToDeg() );
  fillHistogram( kPsi, psi*TMath::RadToDeg() );
  fillHistogram( kt, t );     
}
