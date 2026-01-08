#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/decayAngles.h"

#include "AMPTOOLS_DATAIO/VecPs3PiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

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
  
   bookHistogram( kProd_Ang, new Histogram1D( 50, -PI, PI, "ProdAng", "Production Angle Lab frame [rad.]" ) );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta_GJ", "cos#theta Gottfried-Jackson frame [rad.]" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -PI, PI, "Phi_GJ", "#phi Gottfried-Jackson frame [rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta_HF", "cos#theta Helicity frame [rad.]" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -PI, PI, "Phi_HF", "#phi Helicity frame [rad.]" ) );

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

  // Fixed target
   TLorentzVector target(0,0,0,0.938272);

  
   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector proton = kin->particle( 1 );
   TLorentzVector bach = kin->particle( 2 );
   TLorentzVector vec_daught1 = kin->particle( 3 );
   TLorentzVector vec_daught2 = kin->particle( 4 );
   TLorentzVector piplusL = kin->particle( 5 );

   

   // final meson system P4
   TLorentzVector vec = vec_daught1 + vec_daught2;
   TLorentzVector X = vec + bach;
   TLorentzVector recoil = proton + piplusL;
   
   TLorentzVector proton_ps = proton + bach;
   TLorentzVector recoil_ps = recoil + bach;
   
   
   // We use only default 2-body vector decay (set flag in config file for omega->3pi)
   bool m_3pi = false;  


   // check config file for optional parameters -- we assume here that the first amplitude in the list is a Vec_ps_refl amplitude
   const vector< string > args = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0);
   for(uint ioption=5; ioption<args.size(); ioption++) {
          TString option = args[ioption].c_str();
	  if(option.EqualTo("omega3pi")) m_3pi = true;
   }

   // if the amplitude is a moment, force the 3pi decay as it's currently the only one handled by it
   const vector < string > momentArgs = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0);
   for (uint ioption=0; ioption<momentArgs.size(); ioption++) {
          TString option = momentArgs[ioption].c_str();
      IF(OPTION.EQUALTO("VEC_PS_MOMENT")) M_3PI = TRUE;
   }



   //Momentum transfer
   double Mandt = fabs((target-recoil).M2());
   
   //Calculate production angle in the Gottfried-Jackson frame
   GDouble prod_angle = getPhiProd(beam_polAngle, X, beam, target, 2, true);

   // Calculate decay angles for X in the Gottfried-Jackson frame and for Isobar in the Helicity frame  
   vector <double> thetaPhiAnglesTwoStep = getTwoStepAngles(X, vec, vec_daught1, TLorentzVector(0,0,0,0), beam, target, 2, true);


   GDouble cosTheta = TMath::Cos(thetaPhiAnglesTwoStep[0]);
   GDouble phi = thetaPhiAnglesTwoStep[1];
   GDouble cosThetaH = TMath::Cos(thetaPhiAnglesTwoStep[2]);
   GDouble phiH = thetaPhiAnglesTwoStep[3];
   
   

   //cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, Phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, PhiH );

   fillHistogram( kVecMass, vec.M() );
   fillHistogram( kVecPsMass, X.M() );

   fillHistogram( kt, Mandt );
   fillHistogram( kRecoilMass, recoil.M() );
   fillHistogram( kProtonPsMass, proton_ps.M() );
   fillHistogram( kRecoilPsMass, recoil_ps.M() );

}
