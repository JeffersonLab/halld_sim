#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/vecPsAngles.h"

#include "AMPTOOLS_DATAIO/VecPsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

// Function to check if a string is a valid number
static bool isValidNumber(const string& argInput, double &value){
    char* end = nullptr;
    errno = 0;  // reset global error
    value = std::strtod(argInput.c_str(), &end);

    // Check if 
    // (1) no conversion was performed 
    // (2) there are leftover characters
    // (3) an overflow/underflow occurred   
    if(end == argInput.c_str() || *end != '\0' || errno != 0) {
        return false;  // not a valid number
    }
    // If end points to the end of string, it's fully numeric
    return true;
}

static double parseValidatedNumber(const string& label, const string& argInput){
    double tmpValue = 0.0;
    if(!isValidNumber(argInput, tmpValue)){
      throw std::invalid_argument("Vec_ps_refl: invalid " + label + ": " + argInput);
    }
    return tmpValue;
}

/* Constructor to display FitResults */
VecPsPlotGenerator::VecPsPlotGenerator( const FitResults& results, Option opt ) :
PlotGenerator( results, opt )
{
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
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void
VecPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   // 1st after proton is always the pseudoscalar (bachelor) meson
   TLorentzVector bach = kin->particle( 2 );

   // Default is 2-body vector decay (set flag in config file for omega->3pi)
   bool m_3pi = false;  

   // Min particle index for recoil sum
   int min_recoil = 5; 

   
   // Properly read polarization angle from config file if provided
   double beam_polAngle=0;
   // Check config file for optional parameters -- we assume here that the
   // first amplitude in the list is a Vec_ps_refl amplitude or a Vec_ps_moment
   const vector<string> args =
       cfgInfo()->amplitudeList(reactionName, "", "").at(0)->factors().at(0);
   for(uint ioption = 0; ioption < args.size(); ioption++){
      TString option = args[ioption].c_str();
      if(option.EqualTo("Vec_ps_moment")){
      // Force the 3pi decay as it's currently the only one handled by it
      m_3pi = true;
      break;
     }
     else{
       if(ioption == 6){
         beam_polAngle = parseValidatedNumber("polarization angle", args[6]);
       }
       if(option.EqualTo("omega3pi")){
         m_3pi = true;
       }
     }
   }

   // Compute vector meson from its decay products
   // Make sure the order of daughters is correct in the config file!
   TLorentzVector vec, vec_daught1, vec_daught2;  
   double dalitz_s, dalitz_t, dalitz_u, dalitz_d, dalitz_sc;
   double dalitz_x = 0.0;
   double dalitz_y = 0.0;

   if(m_3pi) {
     // Omega ps proton, omega -> 3pi (6 particles)
     // Omega pi- Delta++, omega -> 3pi (7 particles)
     TLorentzVector pi0 = kin->particle( 3 );
     TLorentzVector pip = kin->particle( 4 );
     TLorentzVector pim = kin->particle( 5 );
     vec = pi0 + pip + pim;
     vec_daught1 = pip;
     vec_daught2 = pim;
	  min_recoil = 6;

	  // Dalitz variables
	  const auto& p2 = pi0;
	  const auto& p3 = pip;
	  const auto& p4 = pim;
     double pdg_m_pi0 = 0.1349766;
     double pdg_m_chargedPi = 0.13957018;
     dalitz_s = (p3 + p4).M2();  // s=M(pip pim)
     dalitz_t = (p2 + p3).M2();  // s=M(pip pi0)
     dalitz_u = (p2 + p4).M2();  // s=M(pim pi0)
     dalitz_d = 2 * (p2 + p3 + p4).M() *
                ((p2 + p3 + p4).M() - ((2 * pdg_m_chargedPi) + pdg_m_pi0));
     dalitz_sc = (1.0 / 3.0) * ((p2 + p3 + p4).M2() +
                                ((2 * (pdg_m_chargedPi * pdg_m_chargedPi)) +
                                 (pdg_m_pi0 * pdg_m_pi0)));
     dalitz_x = sqrt(3.0) * (dalitz_t - dalitz_u) / dalitz_d;
     dalitz_y = 3.0 * (dalitz_sc - dalitz_s) / dalitz_d;

   }
   else {
	  // omega ps proton, omega -> pi0 g (4 particles)
	  // omega pi- Delta++, omega -> pi0 g (5 particles)
	  
	  // (vec 2-body) ps proton, vec 2-body -> pipi, KK (5 particles)
	  // (vec 2-body) pi- Delta++, vec 2-body -> pipi, KK (6 particles)
	  // (vec 2-body) K+ Lambda, vec 2-body -> Kpi (6 particles)
	  vec_daught1 = kin->particle( 3 );
     vec_daught2 = kin->particle( 4 );
     vec = vec_daught1 + vec_daught2;
	  min_recoil = 5;
   }

   // Final meson system P4
   TLorentzVector X = vec + bach;

   TLorentzVector proton_ps = recoil + bach;
   TLorentzVector recoil_ps = proton_ps;
   for(uint i=min_recoil; i<kin->particleList().size(); i++){ 
	   // Add mesons to recoil system (e.g. Delta or Lambda)
	   recoil += kin->particle(i);
	   recoil_ps += kin->particle(i);
   }

   TLorentzVector target(0,0,0,0.938);
   // Helicity coordinate system
   TLorentzVector beamP = beam + target;

   // Calculate decay angles for X in helicity frame (same for all vectors)
   // Change getXDecayAngles to get Gottfried-Jackson angles if needed
   // Note: it also calculates the production angle
   vector <double> xDecayAngles = getXDecayAngles( beam_polAngle, beam, beamP, X, vec);
   // set polarization angle to zero to see shift in Phi_Prod distributions 
   vector <double> xDecayAngles_offset = getXDecayAngles( 0, beam, beamP, X, vec);

   // Calculate vector decay angles (unique for each vector)
   vector <double> vectorDecayAngles;
   if(m_3pi){
    vectorDecayAngles = getVectorDecayAngles( beamP, X, vec,
                                              vec_daught1, vec_daught2);
   }
   else{
    vectorDecayAngles = getVectorDecayAngles( beamP, X, vec,
                                        vec_daught1, TLorentzVector(0,0,0,0));
   }

   double Mandt = fabs((target-recoil).M2());
   double recoil_mass = recoil.M();  

   GDouble cosTheta = TMath::Cos(xDecayAngles[0]);
   GDouble phi = xDecayAngles[1];
   GDouble prod_angle = xDecayAngles[2]; // bigPhi
   GDouble prod_angle_offset = xDecayAngles_offset[2]; // bigPhi not corrected
   GDouble cosThetaH = TMath::Cos(vectorDecayAngles[0]);
   GDouble phiH = vectorDecayAngles[1];
   GDouble lambda = vectorDecayAngles[2];

   //cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kVecPsMass, X.M() );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, phiH );
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kt, Mandt );
   fillHistogram( kRecoilMass, recoil_mass );
   fillHistogram( kProtonPsMass, proton_ps.M() );
   fillHistogram( kRecoilPsMass, recoil_ps.M() );
   fillHistogram( kLambda, lambda );
   fillHistogram( kDalitz, dalitz_x, dalitz_y );
   fillHistogram( kPhi_ProdVsPhi, phi, prod_angle );
   fillHistogram( kPhiOffsetVsPhi, phi, prod_angle_offset );
}
