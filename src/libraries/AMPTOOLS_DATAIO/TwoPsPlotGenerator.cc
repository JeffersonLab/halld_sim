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

  // SDME model
  bookHistogram( kA000,  new Histogram1D( 100, -1.0, 1.0, "A000",  "cos^{2}(#theta)" ) );
  bookHistogram( kA100,  new Histogram1D( 100, -1.0, 1.0, "A100",  "sin(2#theta)cos#phi" ) );
  bookHistogram( kA1m10, new Histogram1D( 100, -1.0, 1.0, "A1m10", "sin^{2}(#theta)cos(2#phi)" ) );

  bookHistogram( kA111,  new Histogram1D( 100, -1.0, 1.0, "A111",  "cos(2#Phi)sin^{2}(#theta)" ) );
  bookHistogram( kA001,  new Histogram1D( 100, -1.0, 1.0, "A001",  "cos(2#Phi)cos^{2}(#theta)" ) );
  bookHistogram( kA101,  new Histogram1D( 100, -1.0, 1.0, "A101",  "cos(2#Phi)sin(2#theta)cos#phi" ) );
  bookHistogram( kA1m11, new Histogram1D( 100, -1.0, 1.0, "A1m11", "cos(2#Phi)sin^{2}(#theta)cos(2#phi)" ) );

  bookHistogram( kA102,  new Histogram1D( 100, -1.0, 1.0, "A102",  "sin(2#Phi)sin(2#theta)sin#phi" ) );
  bookHistogram( kA1m12, new Histogram1D( 100, -1.0, 1.0, "A1m12", "sin(2#Phi)sin^{2}(#theta)sin(2#phi)" ) );

  // K+ polarization model
  bookHistogram( kB1, new Histogram1D( 100, -1.0, 1.0, "B1", "cos#theta_{Y}^{#Lambda(hel)}" ) );
  bookHistogram( kB2, new Histogram1D( 100, -1.0, 1.0, "B2", "cos(2#Phi)" ) );
  bookHistogram( kB3, new Histogram1D( 100, -1.0, 1.0, "B3", "cos(2#Phi)cos#theta_{Y}^{#Lambda(hel)}" ) );
  bookHistogram( kB4, new Histogram1D( 100, -1.0, 1.0, "B4", "sin(2#Phi)cos#theta_{X}^{#Lambda(hel)}" ) );
  bookHistogram( kB5, new Histogram1D( 100, -1.0, 1.0, "B5", "sin(2#Phi)cos#theta_{Z}^{#Lambda(hel)}" ) );

  // Extended Model
  // ---------- X ----------
  bookHistogram( kPx0,    new Histogram1D( 100, -1.0, 1.0, "Px0",    "cos#theta_{X}^{#Lambda(hel)}" ) );
  bookHistogram( kAx000,  new Histogram1D( 100, -1.0, 1.0, "Ax000",  "cos#theta_{X}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAx10c0, new Histogram1D( 100, -1.0, 1.0, "Ax10c0", "cos#theta_{X}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAx10s0, new Histogram1D( 100, -1.0, 1.0, "Ax10s0", "cos#theta_{X}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAx1m1c0,new Histogram1D( 100, -1.0, 1.0, "Ax1m1c0","cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAx1m1s0,new Histogram1D( 100, -1.0, 1.0, "Ax1m1s0","cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPx1,    new Histogram1D( 100, -1.0, 1.0, "Px1",    "cos(2#Phi) cos#theta_{X}^{#Lambda(hel)}" ) );
  bookHistogram( kAx001,  new Histogram1D( 100, -1.0, 1.0, "Ax001",  "cos(2#Phi) cos#theta_{X}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAx10c1, new Histogram1D( 100, -1.0, 1.0, "Ax10c1", "cos(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAx10s1, new Histogram1D( 100, -1.0, 1.0, "Ax10s1", "cos(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAx1m1c1,new Histogram1D( 100, -1.0, 1.0, "Ax1m1c1","cos(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAx1m1s1,new Histogram1D( 100, -1.0, 1.0, "Ax1m1s1","cos(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPx2,    new Histogram1D( 100, -1.0, 1.0, "Px2",    "sin(2#Phi) cos#theta_{X}^{#Lambda(hel)}" ) );
  bookHistogram( kAx002,  new Histogram1D( 100, -1.0, 1.0, "Ax002",  "sin(2#Phi) cos#theta_{X}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAx10c2, new Histogram1D( 100, -1.0, 1.0, "Ax10c2", "sin(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAx10s2, new Histogram1D( 100, -1.0, 1.0, "Ax10s2", "sin(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAx1m1c2,new Histogram1D( 100, -1.0, 1.0, "Ax1m1c2","sin(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAx1m1s2,new Histogram1D( 100, -1.0, 1.0, "Ax1m1s2","sin(2#Phi) cos#theta_{X}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  // ---------- Y ----------
  bookHistogram( kPy0,    new Histogram1D( 100, -1.0, 1.0, "Py0",    "cos#theta_{Y}^{#Lambda(hel)}" ) );
  bookHistogram( kAy000,  new Histogram1D( 100, -1.0, 1.0, "Ay000",  "cos#theta_{Y}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAy10c0, new Histogram1D( 100, -1.0, 1.0, "Ay10c0", "cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAy10s0, new Histogram1D( 100, -1.0, 1.0, "Ay10s0", "cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAy1m1c0,new Histogram1D( 100, -1.0, 1.0, "Ay1m1c0","cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAy1m1s0,new Histogram1D( 100, -1.0, 1.0, "Ay1m1s0","cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPy1,    new Histogram1D( 100, -1.0, 1.0, "Py1",    "cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)}" ) );
  bookHistogram( kAy001,  new Histogram1D( 100, -1.0, 1.0, "Ay001",  "cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAy10c1, new Histogram1D( 100, -1.0, 1.0, "Ay10c1", "cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAy10s1, new Histogram1D( 100, -1.0, 1.0, "Ay10s1", "cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAy1m1c1,new Histogram1D( 100, -1.0, 1.0, "Ay1m1c1","cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAy1m1s1,new Histogram1D( 100, -1.0, 1.0, "Ay1m1s1","cos(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPy2,    new Histogram1D( 100, -1.0, 1.0, "Py2",    "sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)}" ) );
  bookHistogram( kAy002,  new Histogram1D( 100, -1.0, 1.0, "Ay002",  "sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAy10c2, new Histogram1D( 100, -1.0, 1.0, "Ay10c2", "sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAy10s2, new Histogram1D( 100, -1.0, 1.0, "Ay10s2", "sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAy1m1c2,new Histogram1D( 100, -1.0, 1.0, "Ay1m1c2","sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAy1m1s2,new Histogram1D( 100, -1.0, 1.0, "Ay1m1s2","sin(2#Phi) cos#theta_{Y}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  // ---------- Z ----------
  bookHistogram( kPz0,    new Histogram1D( 100, -1.0, 1.0, "Pz0",    "cos#theta_{Z}^{#Lambda(hel)}" ) );
  bookHistogram( kAz000,  new Histogram1D( 100, -1.0, 1.0, "Az000",  "cos#theta_{Z}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAz10c0, new Histogram1D( 100, -1.0, 1.0, "Az10c0", "cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAz10s0, new Histogram1D( 100, -1.0, 1.0, "Az10s0", "cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAz1m1c0,new Histogram1D( 100, -1.0, 1.0, "Az1m1c0","cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAz1m1s0,new Histogram1D( 100, -1.0, 1.0, "Az1m1s0","cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPz1,    new Histogram1D( 100, -1.0, 1.0, "Pz1",    "cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)}" ) );
  bookHistogram( kAz001,  new Histogram1D( 100, -1.0, 1.0, "Az001",  "cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAz10c1, new Histogram1D( 100, -1.0, 1.0, "Az10c1", "cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAz10s1, new Histogram1D( 100, -1.0, 1.0, "Az10s1", "cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAz1m1c1,new Histogram1D( 100, -1.0, 1.0, "Az1m1c1","cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAz1m1s1,new Histogram1D( 100, -1.0, 1.0, "Az1m1s1","cos(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );

  bookHistogram( kPz2,    new Histogram1D( 100, -1.0, 1.0, "Pz2",    "sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)}" ) );
  bookHistogram( kAz002,  new Histogram1D( 100, -1.0, 1.0, "Az002",  "sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)} cos^{2}(#theta)" ) );
  bookHistogram( kAz10c2, new Histogram1D( 100, -1.0, 1.0, "Az10c2", "sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) cos#phi" ) );
  bookHistogram( kAz10s2, new Histogram1D( 100, -1.0, 1.0, "Az10s2", "sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin(2#theta) sin#phi" ) );
  bookHistogram( kAz1m1c2,new Histogram1D( 100, -1.0, 1.0, "Az1m1c2","sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) cos(2#phi)" ) );
  bookHistogram( kAz1m1s2,new Histogram1D( 100, -1.0, 1.0, "Az1m1s2","sin(2#Phi) cos#theta_{Z}^{#Lambda(hel)} sin^{2}(#theta) sin(2#phi)" ) );


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

  // Omega_V
  GDouble cosTheta = angles.CosTheta();
  GDouble Theta = angles.Theta();
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

  // full model
  // 1. Linear Polarization factors
  GDouble g1 = TMath::Cos(2.0 * Phi);
  GDouble g2 = TMath::Sin(2.0 * Phi);

  // 2. Hyperon Decay factors
  GDouble nx = cosThetaX_LambdaHel;
  GDouble ny = cosThetaY_LambdaHel;
  GDouble nz = cosThetaZ_LambdaHel;

  // 3. vector meson simple harmonics
  GDouble b00   = TMath::Cos(Theta) * TMath::Cos(Theta);
  GDouble b11   = TMath::Sin(Theta) * TMath::Sin(Theta);
  GDouble b10c  = TMath::Sin(2.0 * Theta) * TMath::Cos(phi);
  GDouble b10s  = TMath::Sin(2.0 * Theta) * TMath::Sin(phi);
  GDouble b1m1c = TMath::Sin(Theta) * TMath::Sin(Theta) * TMath::Cos(2.0 * phi);
  GDouble b1m1s = TMath::Sin(Theta) * TMath::Sin(Theta) * TMath::Sin(2.0 * phi);


  // SDME model 
  GDouble A000  = b00;   // rho000
  GDouble A100  = b10c;  // rho100
  GDouble A1m10 = b1m1c; // rho1m10
  
  GDouble A111  = g1 * b11;   // rho111
  GDouble A001  = g1 * b00;   // rho001
  GDouble A101  = g1 * b10c;  // rho101
  GDouble A1m11 = g1 * b1m1c; // rho1m11

  GDouble A102  = g2 * b10s;  // rho102
  GDouble A1m12 = g2 * b1m1s; // rho1m12
  
  // K+ Polarization Model
  GDouble B1 = ny; // P
  GDouble B2 = g1; // Sigma
  GDouble B3 = g1 * ny; // alpha * T
  GDouble B4 = g2 * nx; // alpha * Ox
  GDouble B5 = g2 * nz; // alpha * Oz

// Extended Model
  // Y components
  GDouble Py0   =   ny;
  GDouble Ay000 =   ny * b00;
  GDouble Ay10c0 =  ny * b10c;
  GDouble Ay10s0 =  ny * b10s;
  GDouble Ay1m1c0 = ny * b1m1c;
  GDouble Ay1m1s0 = ny * b1m1s;

  GDouble Py1   =   g1 * ny;
  GDouble Ay001 =   g1 * ny * b00;
  GDouble Ay10c1 =  g1 * ny * b10c;
  GDouble Ay10s1 =  g1 * ny * b10s;
  GDouble Ay1m1c1 = g1 * ny * b1m1c;
  GDouble Ay1m1s1 = g1 * ny * b1m1s;

  GDouble Py2   =   g2 * ny;
  GDouble Ay002 =   g2 * ny * b00;
  GDouble Ay10c2 =  g2 * ny * b10c;
  GDouble Ay10s2 =  g2 * ny * b10s;
  GDouble Ay1m1c2 = g2 * ny * b1m1c;
  GDouble Ay1m1s2 = g2 * ny * b1m1s;
  
  // X components
  GDouble Px0   =   nx;
  GDouble Ax000 =   nx * b00;
  GDouble Ax10c0 =  nx * b10c;
  GDouble Ax10s0 =  nx * b10s;
  GDouble Ax1m1c0 = nx * b1m1c;
  GDouble Ax1m1s0 = nx * b1m1s;

  GDouble Px1   =   g1 * nx;
  GDouble Ax001 =   g1 * nx * b00;
  GDouble Ax10c1 =  g1 * nx * b10c;
  GDouble Ax10s1 =  g1 * nx * b10s;
  GDouble Ax1m1c1 = g1 * nx * b1m1c;
  GDouble Ax1m1s1 = g1 * nx * b1m1s;

  GDouble Px2   =   g2 * nx;
  GDouble Ax002 =   g2 * nx * b00;
  GDouble Ax10c2 =  g2 * nx * b10c;
  GDouble Ax10s2 =  g2 * nx * b10s;
  GDouble Ax1m1c2 = g2 * nx * b1m1c;
  GDouble Ax1m1s2 = g2 * nx * b1m1s;

  // Z components
  GDouble Pz0   =   nz;
  GDouble Az000 =   nz * b00;
  GDouble Az10c0 =  nz * b10c;
  GDouble Az10s0 =  nz * b10s;
  GDouble Az1m1c0 = nz * b1m1c;
  GDouble Az1m1s0 = nz * b1m1s;

  GDouble Pz1   =   g1 * nz;
  GDouble Az001 =   g1 * nz * b00;
  GDouble Az10c1 =  g1 * nz * b10c;
  GDouble Az10s1 =  g1 * nz * b10s;
  GDouble Az1m1c1 = g1 * nz * b1m1c;
  GDouble Az1m1s1 = g1 * nz * b1m1s;

  GDouble Pz2   =   g2 * nz;
  GDouble Az002 =   g2 * nz * b00;
  GDouble Az10c2 =  g2 * nz * b10c;
  GDouble Az10s2 =  g2 * nz * b10s;
  GDouble Az1m1c2 = g2 * nz * b1m1c;
  GDouble Az1m1s2 = g2 * nz * b1m1s;

  




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

  fillHistogram( kA000,  A000  );
  fillHistogram( kA100,  A100  );
  fillHistogram( kA1m10, A1m10 );
  fillHistogram( kA111,  A111  );
  fillHistogram( kA001,  A001  );
  fillHistogram( kA101,  A101  );
  fillHistogram( kA1m11, A1m11 );
  fillHistogram( kA102,  A102  );
  fillHistogram( kA1m12, A1m12 );
  fillHistogram( kB1, B1 );
  fillHistogram( kB2, B2 );
  fillHistogram( kB3, B3 );
  fillHistogram( kB4, B4 );
  fillHistogram( kB5, B5 );

  // Extended Model
  fillHistogram( kPx0,    Px0  );
  fillHistogram( kAx000,  Ax000  );
  fillHistogram( kAx10c0, Ax10c0 );
  fillHistogram( kAx10s0, Ax10s0 );
  fillHistogram( kAx1m1c0,Ax1m1c0 );
  fillHistogram( kAx1m1s0,Ax1m1s0 );

  fillHistogram( kPx1,    Px1  );
  fillHistogram( kAx001,  Ax001  );
  fillHistogram( kAx10c1, Ax10c1 );
  fillHistogram( kAx10s1, Ax10s1 );
  fillHistogram( kAx1m1c1,Ax1m1c1 );
  fillHistogram( kAx1m1s1,Ax1m1s1 );

  fillHistogram( kPx2,    Px2  );
  fillHistogram( kAx002,  Ax002  );
  fillHistogram( kAx10c2, Ax10c2 );
  fillHistogram( kAx10s2, Ax10s2 );
  fillHistogram( kAx1m1c2,Ax1m1c2 );
  fillHistogram( kAx1m1s2,Ax1m1s2 );

  fillHistogram( kPy0,    Py0  );
  fillHistogram( kAy000,  Ay000  );
  fillHistogram( kAy10c0, Ay10c0 );
  fillHistogram( kAy10s0, Ay10s0 );
  fillHistogram( kAy1m1c0,Ay1m1c0 );
  fillHistogram( kAy1m1s0,Ay1m1s0 );

  fillHistogram( kPy1,    Py1  );
  fillHistogram( kAy001,  Ay001  );
  fillHistogram( kAy10c1, Ay10c1 );
  fillHistogram( kAy10s1, Ay10s1 );
  fillHistogram( kAy1m1c1,Ay1m1c1 );
  fillHistogram( kAy1m1s1,Ay1m1s1 );

  fillHistogram( kPy2,    Py2  );
  fillHistogram( kAy002,  Ay002  );
  fillHistogram( kAy10c2, Ay10c2 );
  fillHistogram( kAy10s2, Ay10s2 );
  fillHistogram( kAy1m1c2,Ay1m1c2 );
  fillHistogram( kAy1m1s2,Ay1m1s2 );

  fillHistogram( kPz0,    Pz0  );
  fillHistogram( kAz000,  Az000  );
  fillHistogram( kAz10c0, Az10c0 );
  fillHistogram( kAz10s0, Az10s0 );
  fillHistogram( kAz1m1c0,Az1m1c0 );
  fillHistogram( kAz1m1s0,Az1m1s0 );

  fillHistogram( kPz1,    Pz1  );
  fillHistogram( kAz001,  Az001  );
  fillHistogram( kAz10c1, Az10c1 );
  fillHistogram( kAz10s1, Az10s1 );
  fillHistogram( kAz1m1c1,Az1m1c1 );
  fillHistogram( kAz1m1s1,Az1m1s1 );

  fillHistogram( kPz2,    Pz2  );
  fillHistogram( kAz002,  Az002  );
  fillHistogram( kAz10c2, Az10c2 );
  fillHistogram( kAz10s2, Az10s2 );
  fillHistogram( kAz1m1c2,Az1m1c2 );
  fillHistogram( kAz1m1s2,Az1m1s2 );
}
