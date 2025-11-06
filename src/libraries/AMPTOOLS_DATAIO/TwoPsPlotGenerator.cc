#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoPsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/FitResults.h"

TwoPsPlotGenerator::TwoPsPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{

  // // set up the reaction polarization angle map
  // vector< string > reactionVec = reactions();
  // for( auto reac = reactionVec.begin(); reac != reactionVec.end(); ++reac ){

  //   // obtain the polarization angle for this reaction by getting the list of amplitudes
  //   // associated with this reaction -- we know all are SDME amplitudes
  //   // take the 10th argument of the first factor of the first amplitude in the first sum
  //   string ampArgument = { cfgInfo()->amplitudeList( *reac, "", "" ).at(0)->factors().at(0).at(10) };

  //   // pick out the name of the parameter
  //   string parName = ampArgument.substr(1,ampArgument.length()-2);

  //   // here results is of type FitResults which is passed into the constructor of the plot generator
  //   cout << ampArgument.data() << endl;
  //   cout << "parName = " << parName.data() << endl; 
  //   cout << "Angle " << results.parValue( parName ) << endl;
  //   m_reactionAngleMap[*reac] = results.parValue( parName );
  
  // set up the reaction polarization angle map
    vector<string> reactionVec = reactions();
    for (auto& reac : reactionVec) {

        const auto& args = cfgInfo()->amplitudeList(reac, "", "")
                              .at(0)->factors().at(0);

        auto it = std::find_if(args.begin(), args.end(), [](const string& s) {
            string lower = s;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            return lower.find("polangle") != string::npos;
        });

        if (it == args.end()) {
            cerr << "Error: no polAngle argument found for reaction " << reac << endl;
            cerr << "Available arguments: ";
            for (const auto& arg : args) {
                cerr << arg << " ";
            }
            cerr << endl;
            cerr << "Please check your amplitude configuration." << endl;
            std::abort();
        }

        string ampArgument = *it;
        string parName = ampArgument.substr(1, ampArgument.length() - 2);

        cout << ampArgument << endl;
        cout << "parName = " << parName << endl;
        cout << "Angle " << results.parValue(parName) << endl;
        m_reactionAngleMap[reac] = results.parValue(parName);
    }
	createHistograms();
}

TwoPsPlotGenerator::TwoPsPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

TwoPsPlotGenerator::~TwoPsPlotGenerator() {
    // No specific cleanup needed; the base class destructor will handle its cleanup
}


void TwoPsPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( k2PsMass, new Histogram1D( 500, 0.5, 2.0, "M2ps", "Invariant Mass of K #pi") );
  bookHistogram( kLambdaKMass, new Histogram1D( 350, 1, 4.5, "MLambdaK", "Invariant Mass of #Lambda K") );
  bookHistogram( kLambdaPiMass, new Histogram1D( 350, 1, 4.5, "MLambdaPi", "Invariant Mass of #Lambda #pi") );
  bookHistogram( kPiCosTheta, new Histogram1D( 200, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiK,  new Histogram1D( 180, -1*PI, PI, "PhiK",  "#Phi_{K}" ) );
  bookHistogram( kPhiPi, new Histogram1D( 180, -1*PI, PI, "PhiPi", "#Phi_{#pi}" ) );
  bookHistogram( kPhiLambda,  new Histogram1D( 180, -1*PI, PI, "PhiLambda", "#Phi_{#Lambda}" ) );
  
  bookHistogram( kThetaK,  new Histogram1D( 200, 0, 20, "ThetaK",  "#Theta_{K}" ) );
  bookHistogram( kThetaPi, new Histogram1D( 200, 0, 20, "ThetaPi", "#Theta_{#pi}" ) );
  bookHistogram( kThetaLambda, new Histogram1D( 200, 0, 90, "ThetaLambda", "#Theta_{#Lambda}" ) );
  
  bookHistogram( kMomK,  new Histogram1D( 180, 0, 9, "MomK",  "p_{K}" ) );
  bookHistogram( kMomPi, new Histogram1D( 180, 0, 9, "MomPi", "p_{#pi}" ) );
  bookHistogram( kMomLambda, new Histogram1D( 180, 0, 3, "MomLambda", "p_{#Lambda}" ) );
  
  bookHistogram( kPhi_LAB, new Histogram1D( 180, -1*PI, PI, "Phi_LAB", "#Phi Polarization (LAB floor)" ) );
  bookHistogram( kPhi, new Histogram1D( 180, -1*PI, PI, "Phi", "#Phi Polarization (Lab Frame offsetted by PolAngle)" ) );
  bookHistogram( kphi, new Histogram1D( 180, -1*PI, PI, "phi", "#phi (Helicity Frame)" ) );
  bookHistogram( kPsi, new Histogram1D( 180, -1*PI, PI, "psi", "#psi (#phi-#Phi)" ) );
  bookHistogram( kt, new Histogram1D( 500, 0, 3.00, "t", "-t" ) );

  bookHistogram( kLambdaMass, new Histogram1D( 400, 1.05, 1.20, "MLambda", "Invariant Mass of #Lambda" ) );
  bookHistogram( kdaughter1Mass, new Histogram1D( 200, 0.8, 1.2, "Mdaughter1", "Invariant Mass of Proton" ) );
  bookHistogram( kdaughter2Mass, new Histogram1D( 200, 0.0, 0.2, "Mdaughter2", "Invariant Mass of Pi-" ) );


  bookHistogram( kCosThetaX_LambdaHel, new Histogram1D( 200, -1., 1., "cosThetaX_LambdaHel", "cos( #theta_{x} ) in #Lambda helicity frame") );
  bookHistogram( kCosThetaY_LambdaHel, new Histogram1D( 200, -1., 1., "cosThetaY_LambdaHel", "cos( #theta_{y} ) in #Lambda helicity frame") );
  bookHistogram( kCosThetaZ_LambdaHel, new Histogram1D( 200, -1., 1., "cosThetaZ_LambdaHel", "cos( #theta_{z} ) in #Lambda helicity frame") );

  bookHistogram( kCosThetaX_Lambda, new Histogram1D( 200, -1., 1., "cosThetaX_Lambda", "cos( #theta_{x} ) in #Lambda unprimed frame") );
  bookHistogram( kCosThetaY_Lambda, new Histogram1D( 200, -1., 1., "cosThetaY_Lambda", "cos( #theta_{y} ) in #Lambda unprimed frame") );
  bookHistogram( kCosThetaZ_Lambda, new Histogram1D( 200, -1., 1., "cosThetaZ_Lambda", "cos( #theta_{z} ) in #Lambda unprimed frame") );


  bookHistogram( kCosTheta_LambdaHel, new Histogram1D( 200, -1., 1., "cosTheta_LambdaHel", "cos( #theta ) in #Lambda helicity frame") );
  bookHistogram( kPhi_LambdaHel, new Histogram1D( 180, -1*PI, PI, "Phi_LambdaHel", "#Phi in #Lambda helicity frame") );
  bookHistogram( kphi_LambdaHel, new Histogram1D( 180, -1*PI, PI, "phi_LambdaHel", "#phi in #Lambda helicity frame") );
  bookHistogram( kPsi_LambdaHel, new Histogram1D( 180, -1*PI, PI, "psi_LambdaHel", "#psi ( #phi - #Phi ) in #Lambda helicity frame") );
}

void TwoPsPlotGenerator::projectEvent( Kinematics* kin ){

  // this function will make this class backwards-compatible with older versions
  // (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
  // the polariation plane in the lab when multiple orientations are used

  projectEvent( kin, "" );
}

void TwoPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

  double polAngle = m_reactionAngleMap[ reactionName ];
  double protonMass = 0.938272; // GeV/c^2
  double lambdaMass = 1.115683; // GeV/c^2
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector hyperon = kin->particle( 1 ); // Lambda
  TLorentzVector daughter1 = kin->particle( 4 ); // proton
  TLorentzVector daughter2 = kin->particle( 5 ); // pi-
  TLorentzVector p1 = kin->particle( 2 ); // K
  TLorentzVector p2 = kin->particle( 3 ); // Pi

  // res frame
  TLorentzVector recoil = daughter1 + daughter2; // Lambda
  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );
  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  
  // choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();

  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   (p1_res.Vect()).Dot(x),
                     (p1_res.Vect()).Dot(y),
                     (p1_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  
  GDouble phi = angles.Phi();
  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector (1,0,0)
  TVector3 eps_lab(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees 

  // Phi is the azimuthal angle between the photon polarization vector (eps)
  // and the production plane, projected into the transverse (x–y) plane.
  // The second arg of the atan2 ensures the correct sign and quadrant via the right-hand rule (atan2).
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
  GDouble Phi_LAB = atan2(y.Dot(eps_lab), beam.Vect().Unit().Dot(eps_lab.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  // for recoiling particle with different rest mass (hyperon), the formula cannot be simplified
  GDouble t = lambdaMass * lambdaMass - 2 * recoil.E() * protonMass + protonMass * protonMass;

  // Lambda helicity frame
  TLorentzVector target ( 0, 0, 0, 0.9382720813);

  TLorentzRotation LambdaHelBoost( -recoil.BoostVector() );
  TLorentzVector resonance_LambdaHel = LambdaHelBoost * resonance;
  TLorentzVector daughter1_LambdaHel = LambdaHelBoost * daughter1;
  TLorentzVector daughter2_LambdaHel = LambdaHelBoost * daughter2;
  TLorentzVector target_LambdaHel = LambdaHelBoost * target;
  TLorentzVector beam_LambdaHel = LambdaHelBoost * beam;
  
  TVector3 yLambdaHel = (beam_LambdaHel.Vect().Unit().Cross(resonance_LambdaHel.Vect().Unit())).Unit();
  TVector3 zLambdaHel = -1. * resonance_LambdaHel.Vect().Unit();
  TVector3 xLambdaHel = yLambdaHel.Cross(zLambdaHel).Unit();
  TVector3 angles_LambdaHel( (daughter1_LambdaHel.Vect()).Dot(xLambdaHel),
                             (daughter1_LambdaHel.Vect()).Dot(yLambdaHel),
                             (daughter1_LambdaHel.Vect()).Dot(zLambdaHel) );
  GDouble cosTheta_LambdaHel = angles_LambdaHel.CosTheta();
  // project the daughter momentum onto each Lambda helicity axis
  GDouble cosThetaX_LambdaHel = daughter1_LambdaHel.Vect().Unit().Dot( xLambdaHel );
  GDouble cosThetaY_LambdaHel = daughter1_LambdaHel.Vect().Unit().Dot( yLambdaHel );
  GDouble cosThetaZ_LambdaHel = daughter1_LambdaHel.Vect().Unit().Dot( zLambdaHel );
  GDouble phi_LambdaHel = angles_LambdaHel.Phi();
  GDouble psi_LambdaHel = phi_LambdaHel - Phi;
  if(psi_LambdaHel < -1*PI) psi_LambdaHel += 2*PI;
  if(psi_LambdaHel > PI) psi_LambdaHel -= 2*PI;

  
  // traditional "unprimed" frame (Ireland PRL)
  TVector3 zUnprimed = beam_LambdaHel.Vect().Unit();
  TVector3 yUnprimed = yLambdaHel;
  TVector3 xUnprimed = yUnprimed.Cross(zUnprimed).Unit();
  GDouble cosThetaX_Lambda = cos(daughter1_LambdaHel.Vect().Angle(xUnprimed));
  GDouble cosThetaY_Lambda = cos(daughter1_LambdaHel.Vect().Angle(yUnprimed));
  GDouble cosThetaZ_Lambda = cos(daughter1_LambdaHel.Vect().Angle(zUnprimed));

  // Phi angle of the Lambda helivity y axis in both frames
  GDouble Phi_LambdaHel = atan2(yLambdaHel.Dot(eps), beam_LambdaHel.Vect().Unit().Dot(eps.Cross(yLambdaHel)));


  // calls to fillHistogram go here
  
  fillHistogram( k2PsMass, ( resonance ).M() );
  fillHistogram( kLambdaKMass, ( recoil+p1 ).M() );
  fillHistogram( kLambdaPiMass, ( recoil+p2 ).M() );
  fillHistogram( kLambdaMass, ( recoil).M() );
  fillHistogram( kdaughter1Mass, ( daughter1).M() );
  fillHistogram( kdaughter2Mass, ( daughter2).M() );
  fillHistogram( kPiCosTheta, cosTheta );
  fillHistogram( kPhiK,  p1.Phi() );
  fillHistogram( kPhiPi, p2.Phi() );
  fillHistogram( kPhiLambda, recoil.Phi() );
  fillHistogram( kThetaK,  p1.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaPi, p2.Theta()*TMath::RadToDeg() );
  fillHistogram( kThetaLambda, recoil.Theta()*TMath::RadToDeg() );
  fillHistogram( kMomK,  p1.P() );
  fillHistogram( kMomPi, p2.P() );
  fillHistogram( kMomLambda, recoil.P() );
  fillHistogram( kPhi_LAB, Phi_LAB );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );
  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
  fillHistogram( kCosThetaX_LambdaHel, cosThetaX_LambdaHel );
  fillHistogram( kCosThetaY_LambdaHel, cosThetaY_LambdaHel );
  fillHistogram( kCosThetaZ_LambdaHel, cosThetaZ_LambdaHel );
  fillHistogram( kCosTheta_LambdaHel, cosTheta_LambdaHel );
  fillHistogram( kPhi_LambdaHel, Phi_LambdaHel );
  fillHistogram( kphi_LambdaHel, phi_LambdaHel );
  fillHistogram( kPsi_LambdaHel, psi_LambdaHel );
  fillHistogram( kCosThetaX_Lambda, cosThetaX_Lambda );
  fillHistogram( kCosThetaY_Lambda, cosThetaY_Lambda );
  fillHistogram( kCosThetaZ_Lambda, cosThetaZ_Lambda );
}
