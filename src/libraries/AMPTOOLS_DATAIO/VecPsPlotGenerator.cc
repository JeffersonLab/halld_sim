#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_AMPS/omegapiAngles.h"

#include "AMPTOOLS_DATAIO/VecPsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

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
  
   bookHistogram( kVecPsMass, new Histogram1D( 200, 0.6, 2., "MVecPs", "Invariant Mass of Vec+Ps [GeV]") );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#phi [rad.]" ) );
   bookHistogram( kCosThetaH, new Histogram1D( 50, -1., 1., "CosTheta_H", "cos#theta_H" ) );
   bookHistogram( kPhiH, new Histogram1D( 50, -1*PI, PI, "Phi_H", "#phi_H [rad.]" ) );
   bookHistogram( kProd_Ang, new Histogram1D( 50, -1*PI, PI, "Prod_Ang", "Prod_Ang [rad.]" ) );
   bookHistogram( kt, new Histogram1D( 100, 0, 2.0 , "t", "-t" ) );
   bookHistogram( kRecoilMass, new Histogram1D( 100, 0.9, 1.9 , "MRecoil", "Invariant Mass of Recoil [GeV]" ) );
   bookHistogram( kProtonPsMass, new Histogram1D( 100, 0.9, 2.9, "MProtonPs", "Invariant Mass of proton and bachelor Ps [GeV]" ) );
   bookHistogram( kRecoilPsMass, new Histogram1D( 100, 0.9, 2.9, "MRecoilPs", "Invariant Mass of recoil and bachelor Ps [GeV]" ) );
  
}

void
VecPsPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used
  projectEvent( kin, "" );
}

void
VecPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

   //cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector bach = kin->particle( 2 );

   bool m_3pi = true;  // vector decays to 3pi by default (omega)
   int min_recoil = 6; // min particle index for recoil sum

   // check config file for 3pi dalitz parameters -- we assume here that the first amplitude in the list is a Vec_ps_refl amplitude
   int nargs  = cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).size();
   if(nargs == 11) m_3pi = false;

   TLorentzVector vec, vec_daught1, vec_daught2; // compute for each final state below 

   if(m_3pi==1) {
	  // omega ps proton, omega -> 3pi (6 particles)
	  // omega pi- Delta++, omega -> 3pi (7 particles)
          TLorentzVector pi0 = kin->particle( 3 );//omega's pi0
          TLorentzVector pip = kin->particle( 4 );//pi-
          TLorentzVector pim = kin->particle( 5 );//pi+
          vec = pi0 + pip + pim;
          vec_daught1 = pip;
          vec_daught2 = pim;
	  min_recoil = 6;
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

}
