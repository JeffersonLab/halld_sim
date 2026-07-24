#include <algorithm>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/vecPsAngles.h"
#include "AMPTOOLS_DATAIO/VecPsPlotGenerator.h"

#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........ Helper Functions ........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void VecPsPlotGenerator::cacheArgs(){
    
    const auto& reactionList = cfgInfo()->reactionList();
    for( const auto& reaction : reactionList ){
      const std::string reactionName = reaction->reactionName();
      const std::vector<std::string>& fullAmplitudeArgs = 
         cfgInfo()->amplitudeList(reactionName, "", "").at(0)->factors().at(0);

      // Remove the amplitude name since the Vec_ps_refl parser expects only 
      // the arguments of the amplitude
      std::vector<std::string> args(fullAmplitudeArgs.begin() + 1, 
                                    fullAmplitudeArgs.end());
      
      Vec_ps_refl::VecPsReflArgs inputArgs;
      const std::string context = "VecPsPlotGenerator";
      // The first argument in the list is Vec_ps_refl or Vec_ps_moment
      if( fullAmplitudeArgs[0] == "Vec_ps_moment" ){ 
        inputArgs.polAngle  = Vec_ps_refl::parseValidatedNumber( "polAngle", args[1], context );
        inputArgs.omega3pi = true;
        inputArgs.polInfoInPhotonP4 = false;
      }
      else if( fullAmplitudeArgs[0] == "Vec_ps_refl" ){
        inputArgs = Vec_ps_refl::parsedArgs( args , context );
      }
      else{
        throw std::invalid_argument(
          "[ " + context + " ] unrecognized amplitude '" + fullAmplitudeArgs[0] +
          "' (expected 'Vec_ps_refl' or 'Vec_ps_moment')" );
      }
      m_args.emplace(reactionName, inputArgs);
    }
  }

/* Constructor to display FitResults */
VecPsPlotGenerator::VecPsPlotGenerator( const FitResults& results, Option opt ) :
PlotGenerator( results, opt )
{
  cacheArgs();
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
VecPsPlotGenerator::VecPsPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void VecPsPlotGenerator::createHistograms( ) {
  cout << " calls to bookHistogram go here" << endl;

  bookHistogram(kVecPsMass, new Histogram1D(200, 0.6, 2., "MVecPs",
                 ";Invariant Mass of Vec+Ps [GeV]"));
  bookHistogram(kCosTheta, new Histogram1D(50, -1., 1., "CosTheta",
                 ";cos#theta"));
  bookHistogram(kPhi, new Histogram1D(50, -PI, PI, "Phi", ";#phi [rad.]"));
  bookHistogram(kCosThetaH, new Histogram1D(50, -1., 1., "CosTheta_H",
                 ";cos#theta_H"));
  bookHistogram(kPhiH, new Histogram1D(50, -PI, PI, "Phi_H", ";#phi_H [rad.]"));
  bookHistogram(kProd_Ang, new Histogram1D(50, -PI, PI, "Prod_Ang",
                 ";Prod_Ang (#Phi) [rad.]"));
  bookHistogram(kProdOffset, new Histogram1D(50, -PI, PI, "ProdOffset", 
                 ";Prod_Ang (#Phi) Uncorrected [rad.]" ) );
  bookHistogram(kt, new Histogram1D(100, 0, 2.0, "t", ";-t"));
  bookHistogram(kRecoilMass, new Histogram1D(100, 0.9, 1.9, "MRecoil",
                 ";Invariant Mass of Recoil [GeV]"));
  bookHistogram(kProtonPsMass, new Histogram1D(100, 0.9, 2.9, "MProtonPs",
                 ";Invariant Mass of proton and bachelor Ps [GeV]"));
  bookHistogram(kRecoilPsMass, new Histogram1D(100, 0.9, 2.9, "MRecoilPs",
                 ";Invariant Mass of recoil and bachelor Ps [GeV]"));
  bookHistogram(kLambda, new Histogram1D(110, 0.0, 1.1, "Lambda",
                 ";#lambda_{#omega}"));
  bookHistogram(kDalitz, new Histogram2D(100, -2., 2., 100, -2., 2., "Dalitz",
                 ";Dalitz X; Dalitz Y"));
  bookHistogram(kPhi_ProdVsPhi, new Histogram2D(25, -PI, PI, 25, -PI, PI,
                 "Phi_ProdVsPhi", ";#phi [rad.]; Prod_Ang (#Phi) [rad.]" ));
  bookHistogram(kPhiOffsetVsPhi, new Histogram2D(25, -PI, PI, 25, -PI, PI,
                 "PhiOffsetVsPhi",
                 ";#phi [rad.]; Prod_Ang (#Phi) Uncorrected [rad.]" ));
}

void
VecPsPlotGenerator::projectEvent( Kinematics* kin ){

  // This function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polarization plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void
VecPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  const auto& inputArgs = m_args.at(reactionName);
  
  bool polInfoInPhotonP4 = inputArgs.polInfoInPhotonP4;
  double beamPolAngle    = inputArgs.polAngle;
  bool omega3pi          = inputArgs.omega3pi;
  
  //cout << "project event" << endl;
  TLorentzVector readBeam   = kin->particle( 0 );
  TLorentzVector recoil     = kin->particle( 1 );
  // 1st after proton is always the pseudoscalar (bachelor) meson
  TLorentzVector ps         = kin->particle( 2 );

  TVector3 eps;              // beam polarization vector (eps)ilon
  TLorentzVector beam;
  if(polInfoInPhotonP4){
    // When pol info is stored in the photon 4-vector in the tree
    // the energy (pKin[0][0]) and the pz (pKin[0][3]) are used as normal 
    beam.SetPxPyPzE( 0.0, 0.0, readBeam.Pz(), readBeam.E() );
    // while the px and py components store the pol info
    // The values should be stored as px = polFraction*cos(polAngle) 
    // and py = polFraction*sin(polAngle)
    eps.SetXYZ(readBeam.Px(), readBeam.Py(), 0.0); 
    beamPolAngle = eps.Phi();
  }
  else{
    beam = readBeam;
  }

  // Min particle index for recoil sum
  int minRecoil = 5; 

  // Compute vector meson from its decay products
  // Make sure the order of daughters is correct in the config file!
  TLorentzVector vec, vecDaught1, vecDaught2;  
  double dalitzS, dalitzT, dalitzU, dalitzD, dalitzSC;
  double dalitzX = 0.0;
  double dalitzY = 0.0;
  
  if(omega3pi) {
    // Omega ps proton, omega -> 3pi (6 particles):
    // beam proton ps pi0 pip pim
    // Omega pi- Delta++, omega -> 3pi (6 particles):
    // beam delta ps pi0 pip pim
    TLorentzVector pi0 = kin->particle( 3 );
    TLorentzVector pip = kin->particle( 4 );
    TLorentzVector pim = kin->particle( 5 );
    vec = pi0 + pip + pim;
    vecDaught1 = pip;
    vecDaught2 = pim;
    minRecoil = 6;
    
	  // Dalitz variables are included here.
    // They amplitude tends to be multiply in the config file but that's not
    // somthing we know here (can we?)
	  const auto& p2 = pi0;
	  const auto& p3 = pip;
	  const auto& p4 = pim;
    double pdgMassPi0 = 0.1349766;
    double pdgMassPi = 0.13957018;
    dalitzS = (p3 + p4).M2();  // s=M(pip pim)
    dalitzT = (p2 + p3).M2();  // s=M(pip pi0)
    dalitzU = (p2 + p4).M2();  // s=M(pim pi0)
    dalitzD = 2 * (p2 + p3 + p4).M() *
               ((p2 + p3 + p4).M() - ((2 * pdgMassPi) + pdgMassPi0));
    dalitzSC = (1.0 / 3.0) * ((p2 + p3 + p4).M2() +
                               ((2 * (pdgMassPi * pdgMassPi)) +
                                (pdgMassPi0 * pdgMassPi0)));
    dalitzX = sqrt(3.0) * (dalitzT - dalitzU) / dalitzD;
    dalitzY = 3.0 * (dalitzSC - dalitzS) / dalitzD;
  }
   else {
    // Omega ps proton, omega -> gpi0 (5 particles):
    // beam proton ps pi0 gamma
    // Omega pi- Delta++, omega -> gpi (5 particles):
    // beam delta ps pi0 gamma
    // Vec(KK) ps proton, (5 particles)
    // beam proton ps K K
	  // Vec(pipi) pi- Delta++, (5 particles)
    // beam proton ps pi pi
	  // Vec(Kpi) K+ Lambda, (5 particles)
    // beam proton ps K pi
	  vecDaught1 = kin->particle( 3 );
    vecDaught2 = kin->particle( 4 );
    vec = vecDaught1 + vecDaught2;
	  minRecoil = 5;
   }

   // Final meson system P4
   TLorentzVector X = vec + ps;

   TLorentzVector protonPs = recoil + ps;
   TLorentzVector recoilPs = protonPs;
   for(uint i=minRecoil; i<kin->particleList().size(); i++){ 
	   // Add mesons to recoil system (e.g. Delta or Lambda)
	   recoil += kin->particle(i);
	   recoilPs += kin->particle(i);
   }

   TLorentzVector target(0,0,0,0.938);  // proton at rest
   TLorentzVector beamTarget = beam + target;

   // Calculate decay angles for X in helicity frame (same for all vectors)
   // Change getXDecayAngles to get Gottfried-Jackson angles if needed
   // Note: it also calculates the production angle
   vector <double> xDecayAngles = getXDecayAngles( beamPolAngle, beam, beamTarget, X, vec);
   // set polarization angle to zero to see shift in Phi_Prod distributions 
   vector <double> xDecayAnglesOffset = getXDecayAngles( 0, beam, beamTarget, X, vec);

   // Calculate vector decay angles (unique for each vector)
   vector <double> vectorDecayAngles;
   if(omega3pi){
    vectorDecayAngles = getVectorDecayAngles( beamTarget, X, vec,
                                              vecDaught1, vecDaught2);
   }
   else{
    vectorDecayAngles = getVectorDecayAngles( beamTarget, X, vec,
                                        vecDaught1, TLorentzVector(0,0,0,0));
   }

   double momentumTransfer = fabs((target-recoil).M2());
   double recoilMass = recoil.M2();  

   GDouble cosTheta = TMath::Cos(xDecayAngles[0]);
   GDouble phi = xDecayAngles[1];
   GDouble prodAngle = xDecayAngles[2]; // bigPhi
   GDouble prodAngleOffset = xDecayAnglesOffset[2]; // bigPhi not corrected
   GDouble cosThetaH = TMath::Cos(vectorDecayAngles[0]);
   GDouble phiH = vectorDecayAngles[1];
   GDouble lambda = vectorDecayAngles[2];

   //cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kVecPsMass, X.M() );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, phiH );
   fillHistogram( kProd_Ang, prodAngle );
   fillHistogram( kt, momentumTransfer );
   fillHistogram( kRecoilMass, recoilMass );
   fillHistogram( kProtonPsMass, protonPs.M() );
   fillHistogram( kRecoilPsMass, recoilPs.M() );
   fillHistogram( kLambda, lambda );
   fillHistogram( kDalitz, dalitzX, dalitzY );
   fillHistogram( kPhi_ProdVsPhi, phi, prodAngle );
   fillHistogram( kPhiOffsetVsPhi, phi, prodAngleOffset );
}
