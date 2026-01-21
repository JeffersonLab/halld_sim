#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/decayAngles.h"

#include "AMPTOOLS_DATAIO/VecPs3PiPlotGenerator.h"
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
VecPs3PiPlotGenerator::VecPs3PiPlotGenerator( const FitResults& results, Option opt ) :
PlotGenerator( results, opt )
{
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
VecPs3PiPlotGenerator::VecPs3PiPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void VecPs3PiPlotGenerator::createHistograms( ) {
  cout << " calls to bookHistogram go here" << endl;
  
   bookHistogram( kProd_Ang, new Histogram1D( 50, -PI, PI, "ProdAng", "Production Angle [rad.]" ) );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta_GJ", "cos#theta^{[GJ]}" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -PI, PI, "Phi_GJ", "#phi^{[GJ]} [rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta_HF", "cos#theta^{[HF]}" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -PI, PI, "Phi_HF", "#phi^{[HF]} [rad.]" ) );

   bookHistogram( kVecMass, new Histogram1D( 200, 0., 3., "MVec", "m(2#pi)  [GeV]") );
   bookHistogram( kVecPsMass, new Histogram1D( 200, 0.2, 3.2, "MVecPs", "m(3#pi)  [GeV]") );
   bookHistogram( kt, new Histogram1D( 100, 0, 1.0 , "t", "-t  [GeV^{2}]" ) );
   bookHistogram( kRecoilMass, new Histogram1D( 100, 0.9, 1.9 , "ProtonPiplusL_M", "m(p#pi^{+}_{L}) [GeV]" ) );
   
   bookHistogram( kProtonPsMass, new Histogram1D( 200, 0.8, 3.8, "ProtonPiminus_M", "m(p#pi^{-}) [GeV]" ) );
   bookHistogram( kRecoilPsMass, new Histogram1D( 200, 1.0, 4.0, "ProtonPiplusLPiminus_M", "m(p#pi^{+}_{L}#pi^{-}) [GeV]" ) );
 
}

void
VecPs3PiPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void
VecPs3PiPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

   // We work only with a 2-body vector decay 

  
  // Fixed target
   TLorentzVector target(0,0,0,0.938272);

   
   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector proton = kin->particle( 1 );
   TLorentzVector bach = kin->particle( 2 );
   TLorentzVector vec_daught1 = kin->particle( 3 );
   TLorentzVector vec_daught2 = kin->particle( 4 );
   TLorentzVector piplusL = kin->particle( 5 );
   

   // Final state P4 momenta
   TLorentzVector X = vec_daught1 + vec_daught2 + bach;
   TLorentzVector recoil = proton + piplusL;

   //Momenta for the 1st permutation   
   TLorentzVector vec_a = vec_daught1 + vec_daught2;
   TLorentzVector proton_ps_a = proton + bach;
   TLorentzVector recoil_ps_a = recoil + bach;
   
   //Momenta for the 2nd permutation (bach <=> vec_daught1)  
   TLorentzVector vec_b = bach + vec_daught2;
   TLorentzVector proton_ps_b = proton + vec_daught1;
   TLorentzVector recoil_ps_b = recoil + vec_daught1;

   

   // Properly read polarization angle from config file if provided
   double beam_polAngle=0;
   // Check config file for optional parameters -- we assume here that the first amplitude in the list is a Vec_ps_refl amplitude
   const vector< string > args = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0);
   for(uint ioption=5; ioption<args.size(); ioption++) {
          TString option = args[ioption].c_str();
	  if(ioption == 6) beam_polAngle = parseValidatedNumber("polarization angle", args[6]);
   }



   //Momentum transfer
   double Mandt = fabs((target-recoil).M2());
   
   //Calculate production angle in the Gottfried-Jackson frame
   GDouble prod_angle = getPhiProd(beam_polAngle, X, beam, target, 2, true);

   // Calculate decay angles for X in the Gottfried-Jackson frame and for Isobar in the Helicity frame  
   // Angles for the 1st permutation
   vector <double> thetaPhiAnglesTwoStep_a = getTwoStepAngles(X, vec_a, vec_daught1, TLorentzVector(0,0,0,0), beam, target, 2, true);

   // Angles for the 2nd permutation (bach <=> vec_daught1)
   vector <double> thetaPhiAnglesTwoStep_b = getTwoStepAngles(X, vec_b, bach, TLorentzVector(0,0,0,0), beam, target, 2, true);

   
   //Symmetrized angles and masses will be passed as vectors of unit length   
   vector <double> cosTheta_a = {TMath::Cos(thetaPhiAnglesTwoStep_a[0])};
   vector <double> cosTheta_b = {TMath::Cos(thetaPhiAnglesTwoStep_b[0])};

   vector <double> phi_a = {thetaPhiAnglesTwoStep_a[1]};
   vector <double> phi_b = {thetaPhiAnglesTwoStep_b[1]};

   vector <double> cosThetaH_a = {TMath::Cos(thetaPhiAnglesTwoStep_a[2])};
   vector <double> cosThetaH_b = {TMath::Cos(thetaPhiAnglesTwoStep_b[2])};
   
   vector <double> phiH_a = {thetaPhiAnglesTwoStep_a[3]};
   vector <double> phiH_b = {thetaPhiAnglesTwoStep_b[3]};
   
   vector <double> vec_mass_a = {vec_a.M()};
   vector <double> vec_mass_b = {vec_b.M()};
   
   vector <double> protonps_mass_a = {proton_ps_a.M()};
   vector <double> protonps_mass_b = {proton_ps_b.M()};
   
   vector <double> recoilps_mass_a = {recoil_ps_a.M()};
   vector <double> recoilps_mass_b = {recoil_ps_b.M()};

   

   //First, insensitive to the symmetrization quantities
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kt, Mandt );
   fillHistogram( kRecoilMass, recoil.M() );
   fillHistogram( kVecPsMass, X.M() );


   //Now, symmetrized quantities 
   fillHistogram( kCosTheta, cosTheta_a, 0.5 );
   fillHistogram( kCosTheta, cosTheta_b, 0.5 );
 
   fillHistogram( kPhi, phi_a, 0.5 );
   fillHistogram( kPhi, phi_b, 0.5 );
 
   fillHistogram( kCosThetaH, cosThetaH_a, 0.5 );
   fillHistogram( kCosThetaH, cosThetaH_b, 0.5 );
 
   fillHistogram( kPhiH, phiH_a, 0.5 );
   fillHistogram( kPhiH, phiH_b, 0.5 );

   
   fillHistogram( kVecMass, vec_mass_a, 0.5 );
   fillHistogram( kVecMass, vec_mass_b, 0.5 );
  
   fillHistogram( kProtonPsMass, protonps_mass_a, 0.5 );
   fillHistogram( kProtonPsMass, protonps_mass_b, 0.5 );
  
   fillHistogram( kRecoilPsMass, recoilps_mass_a, 0.5 );
   fillHistogram( kRecoilPsMass, recoilps_mass_b, 0.5 );


   
   
}
