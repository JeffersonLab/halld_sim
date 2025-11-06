#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/KHyperonPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

KHyperonPlotGenerator::KHyperonPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{

	vector< string > reactionVec = reactions();
	for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

		// obtain the polarization angle for this reaction by getting the list of amplitudes
		// associated with this reaction -- we know all are SDME amplitudes
		// take the 6th argument of the first factor of the first amplitude in the first sum
		string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(7) };
		
		// pick out the name of the parameter
		string parName = ampArgument.substr(1,ampArgument.length()-2);
		cout<<"angle="<<parName.data()<<" "<<ampArgument.data()<<endl;
	}

	createHistograms();
}

KHyperonPlotGenerator::KHyperonPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void KHyperonPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( kHyperonMass, new Histogram1D( 500, 0.5, 2.0, "MHyp", "Invariant Mass of Hyperon") );

  bookHistogram( kCosThetax, new Histogram1D( 200, -1., 1., "cosThetax", "cos( #theta_{x} ) of Hyperon decay") );
  bookHistogram( kCosThetay, new Histogram1D( 200, -1., 1., "cosThetay", "cos( #theta_{y} ) of Hyperon decay") );
  bookHistogram( kCosThetaz, new Histogram1D( 200, -1., 1., "cosThetaz", "cos( #theta_{z} ) of Hyperon decay") );
  bookHistogram( kCosThetaHyp, new Histogram1D( 200, -1., 1., "cosThetaHyp", "cos( #theta ) of Hyperon decay") );
  bookHistogram( kphiHyp, new Histogram1D( 180, -1*PI, PI, "phiHyp", "#phi of Hyperon decay" ) );
  bookHistogram( kPhiHyp, new Histogram1D( 180, -1*PI, PI, "PhiHyp", "#Phi of Hyperon decay" ) );
}

void KHyperonPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used

  projectEvent( kin, "" );
}

void KHyperonPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){
	
  TLorentzVector beam = kin->particle( 0 ); 
  TLorentzVector k = kin->particle( 1 ); 
  TLorentzVector y1 = kin->particle( 2 ); 
  TLorentzVector y2 = kin->particle( 3 ); 
  
  TLorentzVector hyperon = y1 + y2;
  TLorentzRotation hyperonBoost( -hyperon.BoostVector() );

  fillHistogram( kHyperonMass, hyperon.M() );
	
  TLorentzVector beam_hyperon = hyperonBoost * beam; // beam photon in hyperon rest frame
  TLorentzVector k_hyperon = hyperonBoost * k;       // kaon in hyperon rest frame
  TLorentzVector y1_hyperon = hyperonBoost * y1;     // proton in hyperon rest frame

  // normal to the production plane (formed by beam and kaon)
  TVector3 y = (beam.Vect().Unit().Cross(-k.Vect().Unit())).Unit();

  // choose helicity frame: z-axis opposite kaon in hyperon rest frame
  TVector3 z = -1. * k_hyperon.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles( (y1_hyperon.Vect()).Dot(x),
		   (y1_hyperon.Vect()).Dot(y),
		   (y1_hyperon.Vect()).Dot(z) );

  GDouble cosThetax = cos(y1_hyperon.Vect().Angle(x));
  GDouble cosThetay = cos(y1_hyperon.Vect().Angle(y));
  GDouble cosThetaz = cos(y1_hyperon.Vect().Angle(z));

  fillHistogram( kCosThetax, cosThetax );
  fillHistogram( kCosThetay, cosThetay );
  fillHistogram( kCosThetaz, cosThetaz );

  TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees	
  //double polAngle = m_reactionAngleMap[ reactionName ];
  
  // compute invariant t
  //GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());
  //fillHistogram( kt, -t );      // fill with -t to make positive

  // recoil 
  TLorentzVector target ( 0, 0, 0, 0.9382720813);
  
  TLorentzVector target_hyp = hyperonBoost * target;
  TLorentzVector beam_hyp = hyperonBoost * beam;
  TLorentzVector pion_hyp = hyperonBoost * y1;
  TLorentzVector proton_hyp = hyperonBoost * y2;

  // helicity(?) frame
  TVector3 yHyp = (beam_hyp.Vect().Unit().Cross(-k_hyperon.Vect().Unit())).Unit();  
  TVector3 zHyp = target_hyp.Vect().Unit();
  TVector3 xHyp = yHyp.Cross(zHyp).Unit();
  TVector3 anglesHyp(   (proton_hyp.Vect()).Dot(xHyp),
			(proton_hyp.Vect()).Dot(yHyp),
			(proton_hyp.Vect()).Dot(zHyp) );

  GDouble cosThetaHyp = anglesHyp.CosTheta();  
  GDouble phiHyp = anglesHyp.Phi();
  GDouble PhiHyp = atan2(yHyp.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(yHyp)));
  
  fillHistogram( kCosThetaHyp, cosThetaHyp );
  fillHistogram( kphiHyp, phiHyp );
  fillHistogram( kPhiHyp, PhiHyp );
  
}
