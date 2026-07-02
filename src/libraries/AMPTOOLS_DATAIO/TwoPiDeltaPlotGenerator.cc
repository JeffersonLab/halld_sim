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
    // take the 1th argument of the first factor of the first amplitude in the first sum
    string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(1) };

    // pick out the name of the parameter
    // string parName = ampArgument.substr(1,ampArgument.length()-2);

    // // here results is of type FitResults which is passed into the constructor of the plot generator
    // //cout << "Angle " << results.parValue( parName ) << endl;
    // m_reactionAngleMap[*reac] = results.parValue( parName );
  

    double polAngle = atof(ampArgument.c_str());
    m_reactionAngleMap[*reac] = polAngle;
    cout << "Angle " << polAngle << endl;
  

  }
//   for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

//     string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(10) };
//     string parName = ampArgument.substr(1,ampArgument.length()-2);

//     m_reactionAngleMap[*reac] = results.parValue( parName );
// }


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
  
  bookHistogram( k2PiMass, new Histogram1D( 50, 0.6, 0.95, "M2pi", "M(#pi^{0} #pi^{-}) (GeV)") );
  bookHistogram( kPPipMass, new Histogram1D( 50, 1.1, 1.4, "Mppip", "M(p#pi^{+} )(GeV)") );

  bookHistogram( kPipPimMass, new Histogram1D( 100, 0.2, 2.5, "Mpippim", "M(#pi^{+}#pi^{-} )(GeV)") );
  bookHistogram( kPPi0Mass, new Histogram1D( 100, 1, 4.5, "Mppi0", "M(p#pi^{0}) (GeV)") );
  bookHistogram( kPipPi0Mass, new Histogram1D( 100, 0.1, 2.5, "Mpippi0", "M(#pi^{+}#pi^{0}) (GeV)") );

  bookHistogram( kPPimMass, new Histogram1D( 100, 1, 4.5, "Mppim", "M( p#pi^{-}) (GeV)") );
  bookHistogram( kPi0PimMass, new Histogram1D( 100, 0.1, 1.5, "Mpi0pim", "M(#pi^{0}#pi^{-}) (GeV)") );
  bookHistogram( kPPipPimMass, new Histogram1D( 100, 1, 4.5, "Mppippim", "M(p#pi^{-}#pi^{+})(GeV)") );
  bookHistogram( kPPipPi0Mass, new Histogram1D( 100, 1, 4.5, "Mppippi0", "M( p#pi^{+}#pi^{0}) (GeV)") );
  bookHistogram( kPPimPi0Mass, new Histogram1D( 100, 2.5, 4.5, "Mppimpi0", "M( p#pi^{-}#pi^{0}) (GeV)") );
  bookHistogram( kPipPi0PimMass, new Histogram1D( 100, 0.5, 4.5, "Mpippi0pim", "M( #pi^{-}#pi^{0}#pi^{+}) (GeV)") ); 
  bookHistogram( kPPipPimMass, new Histogram1D( 100, 1, 4.5, "Mppippim", "M( p#pi^{-}#pi^{+}) (GeV)") );
  bookHistogram( kCosTheta_HF_rho, new Histogram1D( 60, -1., 1., "cosTheta_HF_rho", "cos( #theta_{#pi^{-}} ) in HF (restframe #rho)") );
  bookHistogram( kCosTheta_GFJ_rho, new Histogram1D( 60, -1., 1., "cosTheta_GFJ_rho", "cos( #theta_{#pi^{-}} ) in GJ (restframe #rho)") );
  bookHistogram( kCosTheta_GFJ_Delta, new Histogram1D( 60, -1., 1., "cosTheta_GFJ_Delta", "cos( #theta_{#pi^{+}} ) in GJ (restframe #Delta)") );
  bookHistogram( kCosTheta_HF_Delta, new Histogram1D( 60, -1., 1., "cosTheta_HF_Delta", "cos( #theta_{#pi^{+}} ) in HF (restframe #Delta)") );
  bookHistogram( kPhiPiPlus,  new Histogram1D( 180, -180, 180, "PhiPiPlus",  "#Phi_{#pi^{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 180, -180, 180, "PhiPiMinus", "#Phi_{#pi^{-}}" ) );
  bookHistogram( kPhiProton,  new Histogram1D( 180, -180, 180, "PhiProton", "#Phi_{p}" ) );
  bookHistogram( kThetaPiPlus,  new Histogram1D( 200, 0, 20, "ThetaPiPlus",  "#Theta_{#pi^{+}}" ) );
  bookHistogram( kThetaPiMinus, new Histogram1D( 200, 0, 20, "ThetaPiMinus", "#Theta_{#pi^{-}}" ) );
  bookHistogram( kThetaDelta, new Histogram1D( 200, 50, 90, "ThetaDelta", "#Theta_{#Delta}" ) );
  bookHistogram( kMomPiPlus,  new Histogram1D( 60, 0, 0.7, "MomPiPlus",  "p_{#pi^{+}}" ) );
  bookHistogram( kMomPiMinus, new Histogram1D( 60, 0, 9, "MomPiMinus", "p_{#pi^{-}}" ) );
  bookHistogram( kMomPi0, new Histogram1D( 60, 0, 9, "MomPi0", "p_{#pi^{0}}" ) );
  bookHistogram( kMomProton, new Histogram1D( 180, 0, 1, "MomProton", "p_{p}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -180, 180, "Phi", "#Phi" ) );
  bookHistogram( kPhi_GFJ, new Histogram1D( 50, -180, 180, "PhiGFJ", "#Phi_{GJ}" ) );

  bookHistogram( klongMomPiPlus,  new Histogram1D( 60, 0, 0.5, "longMomPiPlus",  "p_{z,#pi^{+}}" ) );
  bookHistogram( klongMomPiMinus, new Histogram1D( 60, 0, 9, "longMomPiMinus", "p_{z,#pi^{-}}" ) );
  bookHistogram( klongMomPi0, new Histogram1D( 60, 0, 9, "longMomPi0", "p_{z,#pi^{0}}" ) );
  bookHistogram( klongMomProton, new Histogram1D( 60, 0, 1.5, "longMomProton", "p_{z,p}" ) );
  
  bookHistogram( kphi_HF_rho, new Histogram1D( 60, -180, 180, "phi_HF_rho", " #phi_{#pi^{-}} in HF (restframe #rho) " ) );
  bookHistogram( kphi_GFJ_rho, new Histogram1D( 60, -180, 180, "phi_GFJ_rho", "#phi_{#pi^{-}} in GJ (restframe #rho)" ) );
  bookHistogram( kphi_HF_Delta, new Histogram1D( 60, -180, 180, "phi_HF_Delta", " #phi_{#pi^{+}} in HF (restframe #Delta)" ) );
  bookHistogram( kphi_GFJ_Delta, new Histogram1D( 60, -180, 180, "phi_GFJ_Delta", "#phi_{#pi^{+}} in GJ (restframe #Delta)" ) );
  bookHistogram( kPsi, new Histogram1D( 20, -180, 180, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 1000, 0, 1.4, "-t", "-t (GeV^2)" ) );
  bookHistogram( kBeamassymetrie, new Histogram1D( 100, 0.16, 0.21 , "Beamass", "#Sigma_{#rho}" ) );
  bookHistogram( kBeamassymetrie_Delta, new Histogram1D( 100, 0.16, 0.21 , "Beamass_Delta", "#Sigma_{#Delta}" ) );
  
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

  TLorentzVector resonance = p2 + p3; 
  TLorentzVector recoil = proton + p1;
  TLorentzVector CM_motion_lab = beam + target;
  TLorentzRotation CMBoost( -CM_motion_lab.BoostVector() );
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p2_res = resonanceBoost * p2;
  TLorentzVector beam_res = resonanceBoost * beam;
  TLorentzVector beam_cm = CMBoost * beam;
  TLorentzVector recoil_cm = CMBoost * recoil;
  TLorentzVector p2_cm = CMBoost * p2;
  TLorentzVector resonance_cm = CMBoost * resonance;

	// logitudinal component



  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  //TVector3 y = (beam_cm.Vect().Unit().Cross(resonance_cm.Vect().Unit())).Unit();   
  //TVector3 y = (beam_res.Vect().Unit().Cross(-recoil_res.Vect().Unit())).Unit();
  
  // choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();

  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   (p2_res.Vect()).Dot(x),
                     (p2_res.Vect()).Dot(y),
                     (p2_res.Vect()).Dot(z) );
 
   // choose GFJ frame: z-axis points in direction of beam photon in restframe

  TVector3 y_GFJ_rho = (beam_res.Vect().Unit().Cross(-recoil_res.Vect().Unit())).Unit();
  
  TVector3 z_GFJ_rho = beam_res.Vect().Unit();

  TVector3 x_GFJ_rho = y_GFJ_rho.Cross(z_GFJ_rho).Unit();
  TVector3 angles_GFJ_rho(   (p2_res.Vect()).Dot(x_GFJ_rho),
                     (p2_res.Vect()).Dot(y_GFJ_rho),
                     (p2_res.Vect()).Dot(z_GFJ_rho) );


  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();

  GDouble cosTheta_GFJ_rho = angles_GFJ_rho.CosTheta();
  GDouble phi_GFJ_rho = angles_GFJ_rho.Phi();
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = -(beam-p3-p2).M2();

  // Calculate angles in GFJ
  TLorentzVector resonance_Delta = proton + p1;   //Delta
  TLorentzVector piminuspi0 = p2 + p3;   // rho 


  
 //Restframe Delta
  TLorentzRotation resonanceBoost_Delta( -resonance_Delta.BoostVector() );

  TLorentzVector target_resDelta = resonanceBoost_Delta * target;
	TLorentzVector beam_resDelta = resonanceBoost_Delta * beam;
	TLorentzVector piminuspi0_resDelta = resonanceBoost_Delta * piminuspi0;
	TLorentzVector p1_resDelta = resonanceBoost_Delta * p1;   // piplus 




 

  TVector3 y_resDelta = (beam_resDelta.Vect().Unit().Cross(piminuspi0_resDelta.Vect().Unit())).Unit();
	TVector3 z_resDelta = target_resDelta.Vect().Unit();

  TVector3 x_resDelta = y_resDelta.Cross(z_resDelta).Unit();
	
	TVector3 angles_piplus_gfj( (p1_resDelta.Vect()).Dot(x_resDelta),
					(p1_resDelta.Vect()).Dot(y_resDelta),
					(p1_resDelta.Vect()).Dot(z_resDelta) );

	GDouble phi_piplus_gfj = angles_piplus_gfj.Phi();
	GDouble cosTheta_piplus_gfj = angles_piplus_gfj.CosTheta();
  GDouble Phi_GFJ = atan2(y_resDelta.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y_resDelta)));  //  _resDelta removed behind beam

  TVector3 z_HF_Delta = -piminuspi0_resDelta.Vect().Unit();
	TVector3 y_HF_Delta = (beam_resDelta.Vect().Cross(piminuspi0_resDelta.Vect())).Unit();
	TVector3 x_HF_Delta = (y_HF_Delta.Cross(z_HF_Delta)).Unit();

  TVector3 angles_HF_Delta( (p1_resDelta.Vect()).Dot(x_HF_Delta),
					(p1_resDelta.Vect()).Dot(y_HF_Delta),
					(p1_resDelta.Vect()).Dot(z_HF_Delta) );

	GDouble phi_HF_Delta= angles_HF_Delta.Phi();
	GDouble cosTheta_HF_Delta = angles_HF_Delta.CosTheta();
  GDouble Beamassymetrie = -0.0881138417353299* 2 * sin(angles_GFJ_rho.Theta())*sin(angles_GFJ_rho.Theta()) -0.231777925342216 *  angles_GFJ_rho.CosTheta() * angles_GFJ_rho.CosTheta();
  GDouble Beamassymetrie_Delta =2 * ((0.0852934984398916 * sin(angles_piplus_gfj.Theta())* sin(angles_piplus_gfj.Theta()))+0.104825277169257*(1/3 +  angles_piplus_gfj.CosTheta() * angles_piplus_gfj.CosTheta()));








 


  // calls to fillHistogram go here
  
  fillHistogram( k2PiMass, ( resonance ).M() );
  fillHistogram( kPPipMass, ( proton+p1 ).M() );
  fillHistogram( kPPimMass, ( proton+p2 ).M() );
  fillHistogram( kPi0PimMass, ( p3+p2 ).M() );
  
  fillHistogram( kCosTheta_HF_rho, cosTheta );
  fillHistogram( kCosTheta_GFJ_Delta, cosTheta_piplus_gfj );
  fillHistogram( kCosTheta_HF_Delta, cosTheta_HF_Delta );
  fillHistogram( kCosTheta_GFJ_rho, cosTheta_GFJ_rho );
  fillHistogram( kPipPimMass,   (p1+p2 ).M());
  fillHistogram( kPPi0Mass,   (proton+p3 ).M());
  fillHistogram( kPipPi0Mass,   (p1+p3 ).M());

  fillHistogram( kPhiPiPlus,  p1.Phi()*TMath::RadToDeg() );
  fillHistogram( kPhiPiMinus, p2.Phi()*TMath::RadToDeg() );
  fillHistogram( kPhiProton, proton.Phi()*TMath::RadToDeg() );
  fillHistogram( kThetaPiPlus,  p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPiMinus, p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaDelta, recoil.Theta()*TMath::RadToDeg() );
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
  fillHistogram( kPhi, Phi*TMath::RadToDeg() );
  fillHistogram( kPhi_GFJ, Phi_GFJ*TMath::RadToDeg() );
  fillHistogram( kphi_HF_rho, phi*TMath::RadToDeg() );
  fillHistogram( kphi_GFJ_Delta, phi_piplus_gfj*TMath::RadToDeg());
  fillHistogram( kphi_HF_Delta, phi_HF_Delta*TMath::RadToDeg());
  fillHistogram( kphi_GFJ_rho, phi_GFJ_rho*TMath::RadToDeg());
  fillHistogram( kPsi, psi*TMath::RadToDeg() );
  fillHistogram( kt, t );      // fill with -t to make positive
  fillHistogram( kBeamassymetrie, Beamassymetrie );
  fillHistogram( kBeamassymetrie_Delta, Beamassymetrie_Delta );
}
