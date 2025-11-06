
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/KStarHyperon.h"


static vector<int> parseIndexList(const string& s) {
  vector<int> indices;
  std::stringstream ss(s);
  int val;
  char comma;
  while (ss >> val) {
      indices.push_back(val);
      if (ss.peek() == ',') ss >> comma;
  }
  return indices;
}

// KStarHyperon::KStarHyperon( const vector<string>& args ) :
//     UserAmplitude<KStarHyperon>( args )
// {
//     assert( args.size() == 10 || args.size() == 11 );

//     // --- NEW: parse K and Lambda indices from args[0], args[1] ---
//     kIndices      = parseIndexList(args[0]);
//     lambdaIndices = parseIndexList(args[1]);
//     assert(kIndices.size() >= 1 && "K must have at least one daughter index");
//     assert(lambdaIndices.size() == 2 && "Lambda must have exactly 2 daughter indices");

//     cout << "K indices: ";
//     for(auto i : kIndices) cout << i << " ";
//     cout << endl;

//     cout << "Lambda indices: ";
//     for(auto i : lambdaIndices) cout << i << " ";
//     cout << endl;

//     // --- Existing amplitude parameters (shifted by +2) ---
//     alpha  = AmpParameter( args[2] );
//     Sigma  = AmpParameter( args[3] );
//     Ox     = AmpParameter( args[4] );
//     P      = AmpParameter( args[5] );
//     T      = AmpParameter( args[6] );
//     Oz     = AmpParameter( args[7] );

//     polAngle = AmpParameter( args[8] );

//     // register parameters
//     registerParameter( alpha );
//     registerParameter( Sigma );
//     registerParameter( Ox );
//     registerParameter( P );
//     registerParameter( T );
//     registerParameter( Oz );
//     registerParameter( polAngle );

//     // polarization handling
//     if(args.size() == 10) {
//         polFraction = atof(args[9].c_str());
//         cout << "Fitting with constant polarization " << polFraction << endl;
//     }
//     else if(args.size() == 11) {
//         polFraction = 0.;
//         TFile* f = new TFile( args[9].c_str() );
//         polFrac_vs_E = (TH1D*)f->Get( args[10].c_str() );
//         assert( polFrac_vs_E != NULL );
//         cout << "Fitting with polarization from " << polFrac_vs_E->GetName() << endl;
//     }
// }

KStarHyperon::KStarHyperon( const vector<string>& args ) :
    UserAmplitude<KStarHyperon>( args )
{
    // now: 2 indices blocks + (alpha, rho111, rho001, Ox, P, T, Oz, polAngle) + (pol spec)
    assert( args.size() == 11 || args.size() == 12 );

    // --- parse K and Lambda indices from args[0], args[1] ---
    kIndices      = parseIndexList(args[0]);
    lambdaIndices = parseIndexList(args[1]);
    assert(kIndices.size() >= 1 && "K must have at least one daughter index");
    assert(lambdaIndices.size() == 2 && "Lambda must have exactly 2 daughter indices");

    cout << "K indices: ";
    for(auto i : kIndices) cout << i << " ";
    cout << endl;

    cout << "Lambda indices: ";
    for(auto i : lambdaIndices) cout << i << " ";
    cout << endl;

    // --- parameters (shifted by +1 because Sigma is now split into two) ---
    alpha   = AmpParameter( args[2] );
    rho111 = AmpParameter( args[3] );   // NEW
    rho001  = AmpParameter( args[4] );   // NEW
    Ox      = AmpParameter( args[5] );
    P       = AmpParameter( args[6] );
    T       = AmpParameter( args[7] );
    Oz      = AmpParameter( args[8] );
    polAngle= AmpParameter( args[9] );

    // register fit parameters 
    registerParameter( alpha );
    registerParameter( rho111 );
    registerParameter( rho001 );
    registerParameter( Ox );
    registerParameter( P );
    registerParameter( T );
    registerParameter( Oz );
    registerParameter( polAngle );

    // (optional) print the current derived initialization for sanity
    cout << "Derived Sigma as (2*rho111 + rho001)" << endl;

    // polarization handling (shifted by +1)
    if(args.size() == 11) {
        // constant polarization in args[10]
        polFraction = atof(args[10].c_str());
        cout << "Fitting with constant polarization " << polFraction << endl;
    }
    else if(args.size() == 12) {
        // histogram-based polarization in args[10] (file) and args[11] (hist name)
        polFraction = 0.;
        TFile* f = new TFile( args[10].c_str() );
        polFrac_vs_E = (TH1D*)f->Get( args[11].c_str() );
        assert( polFrac_vs_E != NULL );
        cout << "Fitting with polarization from " << polFrac_vs_E->GetName() << endl;
    }
}



complex< GDouble >
KStarHyperon::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  GDouble cosThetaX = userVars[kCosThetaX];
  GDouble cosThetaY = userVars[kCosThetaY];
  GDouble cosThetaZ = userVars[kCosThetaZ];
  GDouble phi = polAngle*0.017453293 + userVars[kPhi]; // rotate Phi (in rad)
  GDouble Pgamma = userVars[kPgamma];
  
  // CLAS paper intensity formulation (DOI 10.1103/physrevc.93.065201)
  GDouble I = 1.0 + alpha*cosThetaY*P;
  //I -= Pgamma*cos(2.0*phi) * (Sigma + alpha*cosThetaY*T); 
  I -= Pgamma*cos(2.0*phi) * (2 * rho111 + rho001 + alpha*cosThetaY*T); 
  I += Pgamma*sin(2.0*phi) * alpha * (cosThetaX*Ox + cosThetaZ*Oz);
  
  return complex< GDouble > ( sqrt(fabs(I)) );
}

void
KStarHyperon::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector k;
  for(int idx : kIndices){
      k += TLorentzVector(pKin[idx][1], pKin[idx][2], pKin[idx][3], pKin[idx][0]);
  }
  TLorentzVector y1( 
      pKin[lambdaIndices[0]][1],
      pKin[lambdaIndices[0]][2],
      pKin[lambdaIndices[0]][3],
      pKin[lambdaIndices[0]][0] );
  TLorentzVector y2( 
        pKin[lambdaIndices[1]][1],
        pKin[lambdaIndices[1]][2],
        pKin[lambdaIndices[1]][3],
        pKin[lambdaIndices[1]][0] );
	
  TLorentzVector hyperon = y1 + y2;
  TLorentzRotation hyperonBoost( -hyperon.BoostVector() );
	
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

  userVars[kCosThetaX] = cos(y1_hyperon.Vect().Angle(x));
  userVars[kCosThetaY] = cos(y1_hyperon.Vect().Angle(y));
  userVars[kCosThetaZ] = cos(y1_hyperon.Vect().Angle(z));

  TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees
  userVars[kPhi] = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	
  GDouble Pgamma;
  if(polFraction > 0.) { // for fitting with constant polarization 
    Pgamma = polFraction;
  }
  else{
    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
      Pgamma = 0.;
    }
    else Pgamma = polFrac_vs_E->GetBinContent(bin);
  }
  userVars[kPgamma] = Pgamma;
}

#ifdef GPU_ACCELERATION
void
KStarHyperon::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  // GPUKStarHyperon_exec( dimGrid, dimBlock, GPU_AMP_ARGS, alpha, Sigma, Ox, P, T, Oz,
	// 		   polAngle );
  GPUKStarHyperon_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
          alpha, rho111, rho001, Ox, P, T, Oz, polAngle );
}

#endif //GPU_ACCELERATION
