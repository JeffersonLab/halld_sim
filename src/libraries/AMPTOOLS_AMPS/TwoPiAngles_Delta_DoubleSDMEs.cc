
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_Delta_DoubleSDMEs.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPiAngles_Delta_DoubleSDMEs::TwoPiAngles_Delta_DoubleSDMEs( const vector< string >& args ) :
  UserAmplitude< TwoPiAngles_Delta_DoubleSDMEs >( args )
{

// assert(args.size() == 22 || args.size() == 23 || args.size() == 24 );
frame = string( args[0] );
polAngle = AmpParameter( args[1] );

polFraction = atof(args[2].c_str());
cout << "Fitting with constant polarization " << polFraction << endl;
// r(t)m,m'_lambda,lambda'_alpha      (rt corresponds to rtilde, m corresponds to rho helicity, lambda to Delta helicity, alpha to photon polarization) 
r00_33_0  = AmpParameter( args[3] );

r11_33_0  = AmpParameter( args[4] );
r11_11_0  = AmpParameter( args[5] );

r1m1_33_0 = AmpParameter( args[6] );
r1m1_11_0 = AmpParameter( args[7] );

r10_33_0  = AmpParameter( args[8] );
r10_11_0  = AmpParameter( args[9] );

r00_31_0  = AmpParameter( args[10] );
r00_3m1_0 = AmpParameter( args[11] ); 

r11_31_0  = AmpParameter( args[12] ); 
r11_3m1_0 = AmpParameter( args[13] ); 

r1m1_31_0  = AmpParameter( args[14] ); 
r1m1_3m1_0 = AmpParameter( args[15] );  

r10_31_0  = AmpParameter( args[16] );
r10_3m1_0 = AmpParameter( args[17] );  

rt1m1_31_0  = AmpParameter( args[18] );  
rt1m1_3m1_0 = AmpParameter( args[19] ); 

rt10_31_0   = AmpParameter( args[20] );  
rt10_3m1_0  = AmpParameter( args[21] ); 

//alpha = 1 (same structure)
r00_33_1  = AmpParameter( args[22] );
r00_11_1  = AmpParameter( args[23] );

r11_33_1  = AmpParameter( args[24] );
r11_11_1  = AmpParameter( args[25] );

r1m1_33_1 = AmpParameter( args[26] );
r1m1_11_1 = AmpParameter( args[27] );

r10_33_1  = AmpParameter( args[28] );
r10_11_1  = AmpParameter( args[29] );

r00_31_1  = AmpParameter( args[30] );
r00_3m1_1 = AmpParameter( args[31] );

r11_31_1  = AmpParameter( args[32] );
r11_3m1_1 = AmpParameter( args[33] );

r1m1_31_1  = AmpParameter( args[34] );
r1m1_3m1_1 = AmpParameter( args[35] );

r10_31_1  = AmpParameter( args[36] );
r10_3m1_1 = AmpParameter( args[37] );

rt1m1_31_1  = AmpParameter( args[38] );
rt1m1_3m1_1 = AmpParameter( args[39] );

rt10_31_1   = AmpParameter( args[40] );
rt10_3m1_1  = AmpParameter( args[41] );

//alpha = 2 (different structure)
rt00_31_2  = AmpParameter( args[42] );
rt00_3m1_2 = AmpParameter( args[43] );

rt11_31_2  = AmpParameter( args[44] );
rt11_3m1_2 = AmpParameter( args[45] );

r1m1_33_2 = AmpParameter( args[46] );
r1m1_11_2 = AmpParameter( args[47] );

r1m1_31_2  = AmpParameter( args[48] );
r1m1_3m1_2 = AmpParameter( args[49] );

rt1m1_31_2  = AmpParameter( args[50] );
rt1m1_3m1_2 = AmpParameter( args[51] );

r10_33_2 = AmpParameter( args[52] );
r10_11_2 = AmpParameter( args[53] );

r10_31_2  = AmpParameter( args[54] );
r10_3m1_2 = AmpParameter( args[55] );

rt10_31_2  = AmpParameter( args[56] );
rt10_3m1_2 = AmpParameter( args[57] );

// Parameter r00_11_0 was eliminated due to normalization constraint, so we don't register it);
registerParameter( r00_33_0 );
  
registerParameter( r11_33_0 );
registerParameter( r11_11_0 );

registerParameter( r1m1_33_0 );
registerParameter( r1m1_11_0 );

registerParameter( r10_33_0 );
registerParameter( r10_11_0 );

registerParameter( r00_31_0 );
registerParameter( r00_3m1_0 );

registerParameter( r11_31_0 );
registerParameter( r11_3m1_0 );

registerParameter( r1m1_31_0 );
registerParameter( r1m1_3m1_0 );

registerParameter( r10_31_0 );
registerParameter( r10_3m1_0 );

registerParameter( rt1m1_31_0 );
registerParameter( rt1m1_3m1_0 );

registerParameter( rt10_31_0 );
registerParameter( rt10_3m1_0 );

registerParameter( r00_33_1 );
registerParameter( r00_11_1 );

registerParameter( r11_33_1 );
registerParameter( r11_11_1 );

registerParameter( r1m1_33_1 );
registerParameter( r1m1_11_1 );

registerParameter( r10_33_1 );
registerParameter( r10_11_1 );

registerParameter( r00_31_1 );
registerParameter( r00_3m1_1 );

registerParameter( r11_31_1 );
registerParameter( r11_3m1_1 );

registerParameter( r1m1_31_1 );
registerParameter( r1m1_3m1_1 );

registerParameter( r10_31_1 );
registerParameter( r10_3m1_1 );

registerParameter( rt1m1_31_1 );
registerParameter( rt1m1_3m1_1 );

registerParameter( rt10_31_1 );
registerParameter( rt10_3m1_1 );

//alpha = 2 (different structure)
registerParameter( rt00_31_2 );
registerParameter( rt00_3m1_2 );

registerParameter( rt11_31_2 );
registerParameter( rt11_3m1_2 );

registerParameter( r1m1_33_2 );
registerParameter( r1m1_11_2 );

registerParameter( r1m1_31_2 );
registerParameter( r1m1_3m1_2 );

registerParameter( rt1m1_31_2 );
registerParameter( rt1m1_3m1_2 );

registerParameter( r10_33_2 );
registerParameter( r10_11_2 );

registerParameter( r10_31_2 );
registerParameter( r10_3m1_2 );

registerParameter( rt10_31_2 );
registerParameter( rt10_3m1_2 );

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
TwoPiAngles_Delta_DoubleSDMEs::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {


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


GDouble sinSqTh = sinSqTheta_proton;
GDouble sin2Th = sin2Theta_proton;
GDouble cosSqTh = cosTheta_proton * cosTheta_proton;

GDouble cphi = cos(phi_proton);
GDouble sphi = sin(phi_proton);
GDouble c2phi = cos(2.0 * phi_proton);
GDouble s2phi = sin(2.0 * phi_proton);

GDouble sinSqThPi = sinSqTheta_pim;
GDouble sin2ThPi = sin2Theta_pim;
GDouble cosSqThPi = cosTheta_pim * cosTheta_pim;

GDouble cphiPi = cos(phi_pim);
GDouble sphiPi = sin(phi_pim);
GDouble c2phiPi = cos(2.0 * phi_pim);
GDouble s2phiPi = sin(2.0 * phi_pim);

GDouble r00_11_0 = 0.5 - r00_33_0 - 2.0 * r11_11_0 - 2.0 * r11_33_0;


GDouble W00_0 = r00_33_0 * sinSqTh + r00_11_0 * (1.0/3.0 + cosSqTh);
GDouble Wb00_0 = (2.0 / sqrt(3.0)) * (r00_31_0 * cphi * sin2Th + r00_3m1_0 * c2phi * sinSqTh);

GDouble W11_0 = r11_33_0 * sinSqTh + r11_11_0 * (1.0/3.0 + cosSqTh);
GDouble Wb11_0 = (2.0 / sqrt(3.0)) * (r11_31_0 * cphi * sin2Th + r11_3m1_0 * c2phi * sinSqTh);

GDouble W1m1_0 = r1m1_33_0 * sinSqTh + r1m1_11_0 * (1.0/3.0 + cosSqTh);
GDouble Wb1m1_0 = (2.0 / sqrt(3.0)) * (r1m1_31_0 * cphi * sin2Th + r1m1_3m1_0 * c2phi * sinSqTh);
GDouble Wt1m1_0 = (2.0 / sqrt(3.0)) * (rt1m1_31_0 * sphi * sin2Th + rt1m1_3m1_0 * s2phi * sinSqTh);

GDouble W10_0 = r10_33_0 * sinSqTh + r10_11_0 * (1.0/3.0 + cosSqTh);
GDouble Wb10_0 = (2.0 / sqrt(3.0)) * (r10_31_0 * cphi * sin2Th + r10_3m1_0 * c2phi * sinSqTh);
GDouble Wt10_0 = (2.0 / sqrt(3.0)) * (rt10_31_0 * sphi * sin2Th + rt10_3m1_0 * s2phi * sinSqTh);

GDouble W0 = cosSqThPi * (W00_0 - Wb00_0) + sinSqThPi * (W11_0 - Wb11_0) - sinSqThPi * (c2phiPi * (W1m1_0 - Wb1m1_0) + s2phiPi * Wt1m1_0) - sqrt(2.0) * sin2ThPi * (cphiPi * (W10_0 - Wb10_0) + sphiPi * Wt10_0);
//alpha = 1 (same structure)
GDouble W00_1  = r00_33_1  * sinSqTh + r00_11_1  * (1.0/3.0 + cosSqTh);
GDouble Wb00_1 = (2.0 / sqrt(3.0)) * (r00_31_1  * cphi * sin2Th + r00_3m1_1  * c2phi * sinSqTh);

GDouble W11_1  = r11_33_1  * sinSqTh + r11_11_1  * (1.0/3.0 + cosSqTh);
GDouble Wb11_1 = (2.0 / sqrt(3.0)) * (r11_31_1  * cphi * sin2Th + r11_3m1_1  * c2phi * sinSqTh);

GDouble W1m1_1  = r1m1_33_1  * sinSqTh + r1m1_11_1  * (1.0/3.0 + cosSqTh);
GDouble Wb1m1_1 = (2.0 / sqrt(3.0)) * (r1m1_31_1  * cphi * sin2Th + r1m1_3m1_1  * c2phi * sinSqTh);
GDouble Wt1m1_1 = (2.0 / sqrt(3.0)) * (rt1m1_31_1 * sphi * sin2Th + rt1m1_3m1_1 * s2phi * sinSqTh);

GDouble W10_1  = r10_33_1  * sinSqTh + r10_11_1  * (1.0/3.0 + cosSqTh);
GDouble Wb10_1 = (2.0 / sqrt(3.0)) * (r10_31_1  * cphi * sin2Th + r10_3m1_1  * c2phi * sinSqTh);
GDouble Wt10_1 = (2.0 / sqrt(3.0)) * (rt10_31_1 * sphi * sin2Th + rt10_3m1_1 * s2phi * sinSqTh);

GDouble W1 =cosSqThPi * (W00_1 - Wb00_1) + sinSqThPi * (W11_1 - Wb11_1) - sinSqThPi * (c2phiPi * (W1m1_1 - Wb1m1_1) + s2phiPi * Wt1m1_1) - sqrt(2.0) * sin2ThPi * (cphiPi * (W10_1 - Wb10_1) + sphiPi * Wt10_1);

// alpha = 2 (different structure)
GDouble Wt00_2 = (2.0 / sqrt(3.0)) * (rt00_31_2 * sphi * sin2Th + rt00_3m1_2 * s2phi * sinSqTh);
GDouble Wt11_2 = (2.0 / sqrt(3.0)) * (rt11_31_2 * sphi * sin2Th + rt11_3m1_2 * s2phi * sinSqTh);

GDouble W1m1_2  = r1m1_33_2 * sinSqTh + r1m1_11_2 * (1.0/3.0 + cosSqTh);
GDouble Wb1m1_2 = (2.0 / sqrt(3.0)) * (r1m1_31_2 * cphi * sin2Th + r1m1_3m1_2 * c2phi * sinSqTh);
GDouble Wt1m1_2 = (2.0 / sqrt(3.0)) * (rt1m1_31_2 * sphi * sin2Th + rt1m1_3m1_2 * s2phi * sinSqTh);

GDouble W10_2  = r10_33_2 * sinSqTh + r10_11_2 * (1.0/3.0 + cosSqTh);
GDouble Wb10_2 = (2.0 / sqrt(3.0)) * (r10_31_2 * cphi * sin2Th + r10_3m1_2 * c2phi * sinSqTh);
GDouble Wt10_2 = (2.0 / sqrt(3.0)) * (rt10_31_2 * sphi * sin2Th + rt10_3m1_2 * s2phi * sinSqTh);
GDouble W2 = cosSqThPi * Wt00_2+ sinSqThPi * Wt11_2+ sinSqThPi * (s2phiPi * (W1m1_2 - Wb1m1_2) - c2phiPi * Wt1m1_2) + sqrt(2.0) * sin2ThPi * (sphiPi * (W10_2 - Wb10_2) + cphiPi * Wt10_2); 

GDouble Wpol = W0 - Pgamma * cos(2.0 * bigPhi) * W1 - Pgamma * sin(2.0 * bigPhi) * W2;
return complex< GDouble > ( sqrt(fabs(Wpol)) );
}

void
TwoPiAngles_Delta_DoubleSDMEs::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); // pip
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); // pim
  TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] ); // pi0
  TLorentzVector target ( 0, 0, 0, 0.9382720813);
	// determine boost vectors
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
  TLorentzVector rho_resDelta = DeltaBoost * rho;
  
// Define coordinate systems: 
// y equal in all frames
  TVector3 y = (beam.Vect().Unit().Cross(-Delta.Vect().Unit())).Unit();

  //rho restframe

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
  
  TVector3 angles_proton_GJ_Delta( (proton_resDelta.Vect()).Dot(x_GJ_Delta),
		   (proton_resDelta.Vect()).Dot(y),
		   (proton_resDelta.Vect()).Dot(z_GJ_Delta) );

  
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
  Pgamma = polFraction;
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
TwoPiAngles_Delta_DoubleSDMEs::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUTwoPiAngles_Delta_DoubleSDMEs_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
    r00_33_0, 
    r11_33_0, r11_11_0,
    r1m1_33_0, r1m1_11_0,
    r10_33_0, r10_11_0,
    r00_31_0, r00_3m1_0,
    r11_31_0, r11_3m1_0,
    r1m1_31_0, r1m1_3m1_0,
    r10_31_0, r10_3m1_0,
    rt1m1_31_0, rt1m1_3m1_0,
    rt10_31_0,  rt10_3m1_0,
    r00_33_1, r00_11_1,
    r11_33_1, r11_11_1,
    r1m1_33_1, r1m1_11_1,
    r10_33_1, r10_11_1,
    r00_31_1, r00_3m1_1,
    r11_31_1, r11_3m1_1,
    r1m1_31_1, r1m1_3m1_1,
    r10_31_1, r10_3m1_1,
    rt1m1_31_1, rt1m1_3m1_1,
    rt10_31_1, rt10_3m1_1,
    rt00_31_2, rt00_3m1_2,
    rt11_31_2, rt11_3m1_2,
    r1m1_33_2, r1m1_11_2,
    r1m1_31_2, r1m1_3m1_2,
    rt1m1_31_2, rt1m1_3m1_2,
    r10_33_2, r10_11_2,
    r10_31_2, r10_3m1_2,
    rt10_31_2, rt10_3m1_2
);
}
#endif //GPU_ACCELERATION