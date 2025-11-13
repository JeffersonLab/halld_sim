#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/omegapiAngles.h"

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
  
   bookHistogram( kVecPsMass, new Histogram1D( 200, 0.2, 3.2, "MVecPs", "Invariant Mass of Vec+Ps [GeV]") );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta_I", "cos#theta_I" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi_I", "#phi_I [rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -1*PI, PI, "Phi", "#phi [rad.]" ) );
   bookHistogram( kProd_Ang, new Histogram1D( 50, -1*PI, PI, "Prod_Ang", "Prod_Ang [rad.]" ) );
   bookHistogram( kt, new Histogram1D( 100, 0, 2.0 , "t", "-t" ) );
   bookHistogram( kRecoilMass, new Histogram1D( 100, 0.9, 1.9 , "MRecoil", "Invariant Mass of Recoil [GeV]" ) );
   bookHistogram( kProtonPsMass, new Histogram1D( 200, 0.8, 3.8, "MProtonPs", "Invariant Mass of proton and bachelor pi- [GeV]" ) );
   bookHistogram( kRecoilPsMass, new Histogram1D( 200, 1.0, 4.0, "MRecoilPs", "Invariant Mass of recoil and bachelor pi- [GeV]" ) );
   bookHistogram( kLambda, new Histogram1D( 110, 0.0, 1.1, "Lambda", "#lambda_{#omega}" ) );
   bookHistogram( kDalitz, new Histogram2D( 100, -2., 2., 100, -2., 2., "Dalitz", "Dalitz XY" ) );
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

   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector bach = kin->particle( 3 );

   // default is 2-body vector decay (set flag in config file for omega->3pi)
   bool m_3pi = false;  

   // min particle index for recoil sum
   int min_recoil = 5; 

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
      if(option.EqualTo("Vec_ps_moment")) m_3pi = true;
   }

   TLorentzVector vec, vec_daught1, vec_daught2; // compute for each final state below 
   double dalitz_s, dalitz_t, dalitz_u, dalitz_d, dalitz_sc;
   double dalitzx = 0;
   double dalitzy = 0;

	  // omega ps proton, omega -> pi0 g (4 particles)
	  // omega pi- Delta++, omega -> pi0 g (5 particles)
	  
	  // (vec 2-body) ps proton, vec 2-body -> pipi, KK (5 particles)
	  // (vec 2-body) pi- Delta++, vec 2-body -> pipi, KK (6 particles)
	  // (vec 2-body) K+ Lambda, vec 2-body -> Kpi (6 particles)
	  vec_daught1 = kin->particle( 1 );
          vec_daught2 = kin->particle( 4 );
          vec = vec_daught1 + vec_daught2;
	  min_recoil = 5;
	  dalitzx = 0.;
	  dalitzy = 0.;


   // final meson system P4
   TLorentzVector X = vec + bach;

   TLorentzVector proton_ps = recoil + bach;
   TLorentzVector recoil_ps = proton_ps;
   for(uint i=min_recoil; i<kin->particleList().size(); i++) { 
	// add mesons to recoil system (e.g. Delta or Lambda)
	recoil += kin->particle(i);
	recoil_ps += kin->particle(i);
   }

   // set polarization angle to zero to see shift in Phi_Prod distributions 
   double polAngle = 0;
   TLorentzVector target(0,0,0,0.938);

   // Helicity coordinate system
   TLorentzVector Gammap = beam + target;

   // Calculate decay angles in helicity frame (same for all vectors)
   vector <double> locthetaphi = getomegapiAngles(polAngle, vec, X, beam, Gammap);

   // Calculate vector decay angles (unique for each vector)
   vector <double> locthetaphih;
   if(m_3pi) locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, vec_daught2);
   else locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, TLorentzVector(0,0,0,0));

   double Mandt = fabs((target-recoil).M2());
   double recoil_mass = recoil.M();  

   GDouble cosTheta = TMath::Cos(locthetaphi[0]);
   GDouble Phi = locthetaphi[1];
   GDouble cosThetaH = TMath::Cos(locthetaphih[0]);
   GDouble PhiH = locthetaphih[1];
   GDouble prod_angle = locthetaphi[2];
   GDouble lambda = locthetaphih[2];

   //cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kVecPsMass, X.M() );
   fillHistogram( kCosTheta, cosTheta );
   fillHistogram( kPhi, Phi );
   fillHistogram( kCosThetaH, cosThetaH );
   fillHistogram( kPhiH, PhiH );
   fillHistogram( kProd_Ang, prod_angle );
   fillHistogram( kt, Mandt );
   fillHistogram( kRecoilMass, recoil_mass );
   fillHistogram( kProtonPsMass, proton_ps.M() );
   fillHistogram( kRecoilPsMass, recoil_ps.M() );
   fillHistogram( kLambda, lambda );
   fillHistogram( kDalitz, dalitzx, dalitzy );

}
