
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_only_upperv.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPiAngles_Delta_only_upperv::TwoPiAngles_Delta_only_upperv( const vector< string >& args ) :
  UserAmplitude< TwoPiAngles_Delta_only_upperv >( args )
{
  assert(args.size() == 12 || args.size() == 13 || args.size() == 14 );
	
  rho000  = AmpParameter( args[0] );
  rho100  = AmpParameter( args[1] );
  rho1m10 = AmpParameter( args[2] );
  rho111  = AmpParameter( args[3] );
  rho001  = AmpParameter( args[4] );
  rho101  = AmpParameter( args[5] );
  rho1m11 = AmpParameter( args[6] );
  rho102  = AmpParameter( args[7] );
  rho1m12 = AmpParameter( args[8] );
  frame = string( args[9] );
  polAngle = AmpParameter( args[10] );
  
  
  
  // need to register any free parameters so the framework knows about them
  registerParameter( rho000 );
  registerParameter( rho100 );
  registerParameter( rho1m10 );

  registerParameter( rho111 );
  registerParameter( rho001 );
  // if(args.size() == 12){
  // registerParameter( rho001 );
  // }
  registerParameter( rho101 );
  registerParameter( rho1m11 );

  registerParameter( rho102 );
  registerParameter( rho1m12 );
  
  //new parameter try

  registerParameter( polAngle );

  // }
  // Two possibilities to initialize this amplitude:
  // 1: 11 arguments, fixed polarization
  //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPiAngles <rho000> ... <rho1m12> <polAngle> <polFraction>
  if(args.size() == 13 || args.size() == 12) {
    polFraction = atof(args[11].c_str());
    cout << "Fitting with constant polarization " << polFraction << endl;
  }
  // 2: 12 arguments, read polarization from histogram <hist> in file <rootFile>
  //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPiAngles <rho000> ... <rho1m12> <polAngle> <rootFile> <hist>
//   else if(args.size() == 13) {
//     polFraction = 0.; 
//     TFile* f = new TFile( args[11].c_str() );
//     polFrac_vs_E = (TH1D*)f->Get( args[12].c_str() );
//     assert( polFrac_vs_E != NULL );
//     cout << "Fitting with polarization from " << polFrac_vs_E->GetName() << endl;
//   }
// 
}


complex< GDouble >
TwoPiAngles_Delta_only_upperv::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  GDouble sinSqTheta = userVars[kSinSqTheta];
  GDouble sin2Theta = userVars[kSin2Theta];  
  GDouble cosTheta = userVars[kCosTheta];
  GDouble phi = userVars[kPhi];
  //GDouble bigPhi = polAngle*0.017453293 + userVars[kBigPhi]; // rotate Phi (in rad)
  GDouble bigPhi = userVars[kBigPhi];
  GDouble Pgamma = userVars[kPgamma];
  
 // vector meson production from K. Schilling et. al. last indize = upper

 
 // Use beam assy. to eliminate one parameter

  GDouble W = 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta*cosTheta - sqrt(2.)*rho100*sin2Theta*cos(phi) - rho1m10*sinSqTheta*cos(2.*phi);
	
  W -= Pgamma*cos(2.*bigPhi) * (rho111*sinSqTheta + rho001*cosTheta*cosTheta - sqrt(2.)*rho101*sin2Theta*cos(phi) - rho1m11*sinSqTheta*cos(2.*phi));
	
  W -= Pgamma*sin(2.*bigPhi) * (sqrt(2.)*rho102*sin2Theta*sin(phi) + rho1m12*sinSqTheta*sin(2.*phi));
	
  W *= 3./(4.*PI);

  return complex< GDouble > ( sqrt(fabs(W)) );
}

void
TwoPiAngles_Delta_only_upperv::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] ); 
  TLorentzVector target ( 0, 0, 0, 0.9382720813);
	
  TLorentzVector resonance = p2 + p3;
  TLorentzVector recoil = proton + p1;
  TLorentzVector CM_motion_lab = beam + target;
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );
  TLorentzRotation CMBoost( -CM_motion_lab.BoostVector() );


	
  TLorentzVector beam_res = resonanceBoost * beam;
  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p2_res = resonanceBoost * p2;
  TLorentzVector beam_cm = CMBoost * beam;
  TLorentzVector recoil_cm = CMBoost * recoil;
  TLorentzVector p2_cm = CMBoost * p2;
  TLorentzVector resonance_cm = CMBoost * resonance;

	
  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  //TVector3 y = (beam_cm.Vect().Unit().Cross(resonance_cm.Vect().Unit())).Unit();   
  //TVector3 y = (beam_res.Vect().Unit().Cross(-recoil_res.Vect().Unit())).Unit();

  // choose helicity frame: z-axis opposite recoil Delta in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles( (p2_res.Vect()).Dot(x),
		   (p2_res.Vect()).Dot(y),
		   (p2_res.Vect()).Dot(z) );
  
  //// Calculation angles GFJ frame


	
  // normal to the production plane
  TVector3 y_GFJ = (beam_res.Vect().Unit().Cross(-recoil_res.Vect().Unit())).Unit();

  // choose GFJ frame: z-axis points in direction of beam vector in rho rest frame
  TVector3 z_GFJ = beam_res.Vect().Unit();
  TVector3 x_GFJ = y_GFJ.Cross(z_GFJ).Unit();
  TVector3 angles_GFJ( (p2_res.Vect()).Dot(x_GFJ),
	(p2_res.Vect()).Dot(y_GFJ),
	(p2_res.Vect()).Dot(z_GFJ) );
  
  if(frame=="Hel"){
  userVars[kCosTheta]   = angles.CosTheta();
  userVars[kSinSqTheta] = sin(angles.Theta())*sin(angles.Theta());
  userVars[kSin2Theta]  = sin(2.*angles.Theta());
  userVars[kPhi] = angles.Phi();
  }
  else{
  userVars[kCosTheta]   = angles_GFJ.CosTheta();
  userVars[kSinSqTheta] = sin(angles_GFJ.Theta())*sin(angles_GFJ.Theta());
  userVars[kSin2Theta]  = sin(2.*angles_GFJ.Theta());
  userVars[kPhi] = angles_GFJ.Phi();
  }  
  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); 
  userVars[kBigPhi] = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	
  // vector meson production from K. Schilling et. al.
  GDouble Pgamma;
  if(polFraction > 0.) { // for fitting with constant polarization 
    Pgamma = polFraction;
  }
  // else{
  //   int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
  //   if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
  //     Pgamma = 0.;
  //   }
  //   else Pgamma = polFrac_vs_E->GetBinContent(bin);
  // }
  userVars[kPgamma] = Pgamma;
}

#ifdef GPU_ACCELERATION
void
TwoPiAngles_Delta_only_upperv::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUTwoPiAngles_Delta_only_upperv_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
			   rho000, rho100, rho1m10,
			   rho111, rho001, rho101,
			   rho1m11, rho102, rho1m12,
			   polAngle);
}

#endif //GPU_ACCELERATION