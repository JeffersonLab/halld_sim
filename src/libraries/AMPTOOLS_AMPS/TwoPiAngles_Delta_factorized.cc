#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_factorized.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPiAngles_Delta_factorized::TwoPiAngles_Delta_factorized( const vector< string >& args ) :
  UserAmplitude< TwoPiAngles_Delta_factorized >( args )
{

// assert(args.size() == 22 || args.size() == 23 || args.size() == 24 );
frame = string( args[0] );
polAngle = AmpParameter( args[1] );

polFraction = atof(args[2].c_str());
cout << "Fitting with constant polarization " << polFraction << endl;

rho000  = AmpParameter( args[3] );
rho100  = AmpParameter( args[4] );
rho1m10 = AmpParameter( args[5] );
rho001 = AmpParameter( args[6] );
rho111  = AmpParameter( args[7] );
rho101  = AmpParameter( args[8] );
rho1m11 = AmpParameter( args[9] );
rho102  = AmpParameter( args[10] );
rho1m12 = AmpParameter( args[11] );

delta_rho011  = AmpParameter( args[12] );
delta_rho031  = AmpParameter( args[13] );
delta_rho03m1 = AmpParameter( args[14] );

registerParameter( rho000 );
registerParameter( rho100 );
registerParameter( rho1m10 );
registerParameter( rho001 );
registerParameter( rho111 );
registerParameter( rho101 );
registerParameter( rho1m11 );

registerParameter( rho102 );
registerParameter( rho1m12 );
registerParameter( delta_rho011 );
registerParameter( delta_rho031 );
registerParameter( delta_rho03m1 );

registerParameter( polAngle );

  // }
  // Two possibilities to initialize this amplitude:
  // 1: 11 arguments, fixed polarization
  //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPiAngles <rho000> ... <rho1m12> <polAngle> <polFraction>
  
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
TwoPiAngles_Delta_factorized::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {


//GDouble bigPhi = polAngle*0.017453293 + userVars[kBigPhi]; // rotate Phi (in rad)
GDouble bigPhi = userVars[kBigPhi];
GDouble Pgamma = userVars[kPgamma];
GDouble cosTheta_pim = userVars[kCosTheta_pim];
GDouble sinSqTheta_pim = userVars[kSinSqTheta_pim];
GDouble sin2Theta_pim = userVars[kSin2Theta_pim];
GDouble phi_pim = userVars[kPhi_pim];
GDouble cosTheta_proton = userVars[kCosTheta_proton];
GDouble sinSqTheta_proton = userVars[kSinSqTheta_proton];
GDouble sin2Theta_proton = userVars[kSin2Theta_proton];
GDouble phi_proton = userVars[kPhi_proton];

  
 // vector meson production from K. Schilling et. al. last indize = upper



// GDouble rho001_val;
// rho001_val = 2*(delta_rho111 + delta_rho133 - rho111); 
// Use beam assy. to eliminate one parameter

GDouble W0_meson = 3.*(0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta_pim*cosTheta_pim - sqrt(2.)*rho100*sin2Theta_pim*cos(phi_pim) - rho1m10*sinSqTheta_pim*cos(2.*phi_pim)); 
GDouble W1_meson= 3.*(rho111*sinSqTheta_pim + rho001*cosTheta_pim*cosTheta_pim - sqrt(2.)*rho101*sin2Theta_pim*cos(phi_pim) - rho1m11*sinSqTheta_pim*cos(2.*phi_pim)); 
GDouble W2_meson = 3. * (sqrt(2.)*rho102*sin2Theta_pim*sin(phi_pim) + rho1m12*sinSqTheta_pim*sin(2.*phi_pim)); 

GDouble W0_baryon = 3.*(0.5 - delta_rho011)*sinSqTheta_proton + delta_rho011*(1.+3.*cosTheta_proton*cosTheta_proton) - 2.*TMath::Sqrt(3.)*delta_rho031*cos(phi_proton)*sin2Theta_proton - 2.*TMath::Sqrt(3.)*delta_rho03m1*cos(2.*phi_proton)*sinSqTheta_proton;

GDouble Wpol = W0_meson*W0_baryon - Pgamma * cos(2.0 * bigPhi) *( W1_meson*W0_baryon) - Pgamma * sin(2.0 * bigPhi) *( W2_meson*W0_baryon);
Wpol *= 1/(4.*PI);



return complex< GDouble > ( sqrt(fabs(Wpol)) );
}

void
TwoPiAngles_Delta_factorized::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); // pip
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); // pim
  TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] ); // pi0
  TLorentzVector target ( 0, 0, 0, 0.9382720813);
	// determine boost vectors
  TLorentzVector recoil = proton + p1;
  TLorentzVector rho = p2 + p3;
  TLorentzVector Delta = proton + p1;
  TLorentzVector CM_motion_lab = beam + target;
  TLorentzVector piminus = p2;

  TLorentzRotation rhoBoost( -rho.BoostVector() );
  TLorentzRotation CMBoost( -CM_motion_lab.BoostVector() );
  TLorentzRotation DeltaBoost( -Delta.BoostVector() );
  
//bost in rho restframe
  TLorentzVector beam_resrho = rhoBoost * beam;
  TLorentzVector Delta_resrho = rhoBoost * Delta;
  TLorentzVector piminus_resrho = rhoBoost * piminus;
  TLorentzVector proton_resrho = rhoBoost * proton;
//boost in CM frame
  TLorentzVector beam_cm = CMBoost * beam;
  TLorentzVector Delta_cm = CMBoost * Delta;
  TLorentzVector piminus_cm = CMBoost * piminus;
  TLorentzVector rho_cm = CMBoost * rho;
  
//boost in recoil frame
  TLorentzVector target_resDelta = DeltaBoost * target;
  TLorentzVector proton_resDelta = DeltaBoost * proton;
  TLorentzVector piplus_resDelta = DeltaBoost * p1;
  TLorentzVector rho_resDelta = DeltaBoost * rho;
  
  // Define coordinate systems: 
  // y equal in all frames
  TVector3 y = (beam.Vect().Unit().Cross(-Delta.Vect().Unit())).Unit();

  //rho system
  
  //Helicity frame
  TVector3 z_hel_rho = -1. * Delta_resrho.Vect().Unit();
  TVector3 x_hel_rho = y.Cross(z_hel_rho).Unit();

  // GJ frame
  TVector3 z_GJ_rho = beam_resrho.Vect().Unit();
  TVector3 x_GJ_rho = y.Cross(z_GJ_rho).Unit();

  //Delta system

  // Helicity frame
  TVector3 z_hel_Delta = -rho_resDelta.Vect().Unit();
	TVector3 x_hel_Delta = (y.Cross(z_hel_Delta)).Unit();

  // GJ frame
  TVector3 z_GJ_Delta = target_resDelta.Vect().Unit();
  TVector3 x_GJ_Delta = y.Cross(z_GJ_Delta).Unit();

  //Angle calculation
  TVector3 angles_piminus_hel_rho( (piminus_resrho.Vect()).Dot(x_hel_rho),
		   (piminus_resrho.Vect()).Dot(y),
		   (piminus_resrho.Vect()).Dot(z_hel_rho) );
  
  TVector3 angles_piminus_GJ_rho( (piminus_resrho.Vect()).Dot(x_GJ_rho),
		   (piminus_resrho.Vect()).Dot(y),
		   (piminus_resrho.Vect()).Dot(z_GJ_rho) );

  

  TVector3 angles_proton_hel_Delta( (proton_resDelta.Vect()).Dot(x_hel_Delta),
		   (proton_resDelta.Vect()).Dot(y),
		   (proton_resDelta.Vect()).Dot(z_hel_Delta) );
  TVector3 angles_piplus_hel_Delta( (piplus_resDelta.Vect()).Dot(x_hel_Delta),
		   (piplus_resDelta.Vect()).Dot(y),
		   (piplus_resDelta.Vect()).Dot(z_hel_Delta) );
  
  TVector3 angles_proton_GJ_Delta( (proton_resDelta.Vect()).Dot(x_GJ_Delta),
		   (proton_resDelta.Vect()).Dot(y),
		   (proton_resDelta.Vect()).Dot(z_GJ_Delta) );

  TVector3 angles_piplus_GJ_Delta( (piplus_resDelta.Vect()).Dot(x_GJ_Delta),
  (piplus_resDelta.Vect()).Dot(y),
  (piplus_resDelta.Vect()).Dot(z_GJ_Delta) );
  
  if(frame=="Hel"){

  userVars[kCosTheta_pim]   = angles_piminus_hel_rho.CosTheta();
  userVars[kSinSqTheta_pim] = sin(angles_piminus_hel_rho.Theta())*sin(angles_piminus_hel_rho.Theta());
  userVars[kSin2Theta_pim]  = sin(2.*angles_piminus_hel_rho.Theta());
  userVars[kPhi_pim] = angles_piminus_hel_rho.Phi();

  userVars[kCosTheta_proton]   = angles_proton_hel_Delta.CosTheta();
  userVars[kSinSqTheta_proton] = sin(angles_proton_hel_Delta.Theta())*sin(angles_proton_hel_Delta.Theta());
  userVars[kSin2Theta_proton]  = sin(2.*angles_proton_hel_Delta.Theta());
  userVars[kPhi_proton] = angles_proton_hel_Delta.Phi();
  }
  else{
  userVars[kCosTheta_pim]   = angles_piminus_GJ_rho.CosTheta();
  userVars[kSinSqTheta_pim] = sin(angles_piminus_GJ_rho.Theta())*sin(angles_piminus_GJ_rho.Theta());
  userVars[kSin2Theta_pim]  = sin(2.*angles_piminus_GJ_rho.Theta());
  userVars[kPhi_pim] = angles_piminus_GJ_rho.Phi();

  userVars[kCosTheta_proton]   = angles_proton_GJ_Delta.CosTheta();
  userVars[kSinSqTheta_proton] = sin(angles_proton_GJ_Delta.Theta())*sin(angles_proton_GJ_Delta.Theta());
  userVars[kSin2Theta_proton]  = sin(2.*angles_proton_GJ_Delta.Theta());
  userVars[kPhi_proton] = angles_proton_GJ_Delta.Phi();
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
TwoPiAngles_Delta_factorized::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUTwoPiAngles_Delta_factorized_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
    rho000,
    rho100,
    rho1m10,
    rho001,
    rho111,
    rho101,
    rho1m11,
    rho102,
    rho1m12,
    delta_rho011,
    delta_rho031,
    delta_rho03m1
);
}
#endif //GPU_ACCELERATION