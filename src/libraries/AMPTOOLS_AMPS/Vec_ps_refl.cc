
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
//#include "AMPTOOLS_AMPS/omegapiAngles.h"
#include "AMPTOOLS_AMPS/decayAngles.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"

#include "UTILITIES/BeamProperties.h"

Vec_ps_refl::Vec_ps_refl( const vector< string >& args ) :
UserAmplitude< Vec_ps_refl >( args )
{
  //assert( args.size() == 11 );
  
  
  m_j = atoi( args[0].c_str() ); // resonance spin J
  m_m = atoi( args[1].c_str() ); // spin projection (Lambda)
  m_l = atoi( args[2].c_str() ); // partial wave L
  m_r = atoi( args[3].c_str() ); // real (+1) or imaginary (-1)
  m_s = atoi( args[4].c_str() ); // sign for polarization in amplitude

  // default polarization information stored in tree
  m_polInTree = true;

  // default is 2-body vector decay (set flag in config file for omega->3pi)
  m_3pi = false; 

  // 5 possibilities to initialize this amplitude:
  // (with <J>: total spin, <m>: spin projection, <l>: partial wave, <r>: +1/-1 for real/imaginary part; <s>: +1/-1 sign in P_gamma term)

  // loop over any additional amplitude arguments to change defaults
  for(uint ioption=5; ioption<args.size(); ioption++) {
	  TString option = args[ioption].c_str();

	  // polarization provided in configuration file
	  if(ioption==5 && option.IsFloat()) {
		  m_polInTree = false;
		  polAngle = atof(args[5].c_str());
	  
		  TString polOption = args[6].c_str();
		  if(polOption.IsFloat()) polFraction = atof(polOption.Data());
		  else if(polOption.Contains(".root")) {
			  polFraction = 0.;
			  TFile* f = new TFile( polOption );
			  polFrac_vs_E = (TH1D*)f->Get( args[7].c_str() );
			  assert( polFrac_vs_E != NULL );
		  }
		  else {
			  cout << "ERROR: Vec_ps_refl beam polarization not set" <<endl;
			  assert(0);
		  }
	  }

	  // other options should be strings
	  if(option.EqualTo("omega3pi")) m_3pi = true;

  }

  // make sure values are reasonable
  assert( abs( m_m ) <= m_j );
  // m_r = +1 for real
  // m_r = -1 for imag
  assert( abs( m_r ) == 1 );
  // m_s = +1 for 1 + Pgamma
  // m_s = -1 for 1 - Pgamma
  assert( abs( m_s ) == 1 );
  
}

void
Vec_ps_refl::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector beam;
  TVector3 eps;
  double beam_polFraction;
  double beam_polAngle;

  if(m_polInTree){
    beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // beam polarization vector;

    beam_polFraction = eps.Mag();
    beam_polAngle = eps.Phi()*TMath::RadToDeg();
  }
  else {
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    beam_polAngle = polAngle;
    
    if(polFraction > 0.) { // for fitting with fixed polarization
	    beam_polFraction = polFraction;
    }
    else { // for fitting with polarization vs E_gamma from input histogram 
	    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
	    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		    beam_polFraction = 0.;
	    } else 
		    beam_polFraction = polFrac_vs_E->GetBinContent(bin);
    }
  }
  
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 

  // common vector and pseudoscalar P4s
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); // 1st after proton
  TLorentzVector vec, vec_daught1, vec_daught2; // compute for each final state below 

  // omega ps proton, omega -> 3pi (6 particles)
  // omega pi- Delta++, omega -> 3pi (7 particles)
  if(m_3pi) {
	  TLorentzVector pi0(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  TLorentzVector pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  TLorentzVector pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
	  vec = pi0 + pip + pim;
	  vec_daught1 = pip;
	  vec_daught2 = pim;
  }
  else {
	  // omega ps proton, omega -> pi0 g (4 particles)
	  // omega pi- Delta++, omega -> pi0 g (5 particles)
	  
	  // (vec 2-body) ps proton, vec 2-body -> pipi, KK (5 particles)
	  // (vec 2-body) pi- Delta++, vec 2-body -> pipi, KK (6 particles)
	  // (vec 2-body) K+ Lambda, vec 2-body -> Kpi (6 particles)
	  vec_daught1 = TLorentzVector(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  vec_daught2 = TLorentzVector(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  vec = vec_daught1 + vec_daught2;
  }

  // final meson system P4
  TLorentzVector X = vec + ps;

  //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  TLorentzVector target(0,0,0,0.938);
  
  TLorentzVector Gammap = beam + target;

  // Calculate decay angles in helicity or Gottfried-Jackson frame
  // set beam polarization angle to 0 degrees; apply diamond orientation in calcAmplitude
  double locPhiProd = getPhiProd( 0, X, beam, target, 2, true );
  //vector <double> locthetaphi = getomegapiAngles(0, vec, X, beam, Gammap);
  vector< double > locDecayAngles;
  if(m_3pi) locDecayAngles = getTwoStepAngles( X, vec, vec_daught1, vec_daught2, beam, target, 2, true ); 
  else locDecayAngles = getTwoStepAngles( X, vec, vec_daught1, TLorentzVector(0,0,0,0), beam, target, 2, true ); 


  // Calculate vector decay angles (unique for each vector)
  //vector <double> locthetaphih;
  //if(m_3pi) locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, vec_daught2);
  //else locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, TLorentzVector(0,0,0,0));

  userVars[uv_cosTheta] = TMath::Cos(locDecayAngles[0]);
  userVars[uv_Phi] = locDecayAngles[1];

  userVars[uv_cosThetaH] = TMath::Cos(locDecayAngles[2]);
  userVars[uv_PhiH] = locDecayAngles[3];

  userVars[uv_prod_Phi] = locPhiProd;

  userVars[uv_MX] = X.M();
  userVars[uv_MVec] = vec.M();
  userVars[uv_MPs] = ps.M();

  userVars[uv_beam_polFraction] = beam_polFraction;
  userVars[uv_beam_polAngle] = beam_polAngle;
/*  
  // Use these to test angular calculation for omegapiDelta
  TLorentzVector pi0 = vec - vec_daught1 - vec_daught2;
  TLorentzVector pip = vec_daught1;
  TLorentzVector pim = vec_daught2;
  TLorentzVector pi_recoil( pKin[6][1], pKin[6][2], pKin[6][3], pKin[6][0] );
  double s_mand = Gammap * Gammap;
  TLorentzVector sqrt_t = X - beam;
  TLorentzVector deltaRecoil = recoil + pi_recoil;

  cout << "Four-momenta in lab frame: (Px, Py, Pz, E)" << endl;
  cout << "beam photon: (" << beam[0] << ", " << beam[1] << ", " << beam[2] << ", " << beam[3] << ")" << endl;
  cout << "recoil proton: (" << recoil[0] << ", " << recoil[1] << ", " << recoil[2] << ", " << recoil[3] << ")" << endl;
  cout << "bachelor pion: (" << ps[0] << ", " << ps[1] << ", " << ps[2] << ", " << ps[3] << ")" << endl;
  cout << "pi0 from omega: (" << pi0[0] << ", " << pi0[1] << ", " << pi0[2] << ", " << pi0[3] << ")" << endl;
  cout << "pi+ from omega: (" << pip[0] << ", " << pip[1] << ", " << pip[2] << ", " << pip[3] << ")" << endl;
  cout << "pi- from omega: (" << pim[0] << ", " << pim[1] << ", " << pim[2] << ", " << pim[3] << ")" << endl;
  cout << "recoil pion: (" << pKin[6][1] << ", " << pKin[6][2] << ", " << pKin[6][3] << ", " << pKin[6][0] << ")" << endl;
  cout << "Mandelstam variables:" << endl;
  cout << "s = " << s_mand << endl;
  cout << "t = " << sqrt_t*sqrt_t << endl;
  cout << "m2_omega = " << vec.M2() << endl;
  cout << "m2_X = " << X.M2() << endl;
  cout << "m2_pi1pi2 = " << (pi0 + pip).M2() << endl;
  cout << "m2_pi2pi3 = " << (pip + pim).M2() << endl;
  cout << "m2_ppip = " << ( pi_recoil + recoil ).M2() << endl;
  cout << "angles:" << endl;
  cout << "cos(theta_omega) = " << userVars[uv_cosTheta] << ", phi_omega = " << userVars[uv_Phi] << endl;
  cout << "cos(theta_H) = " << userVars[uv_cosThetaH] << ", phi_H = " << userVars[uv_PhiH] << endl;
  cout << "cos(theta_omega)_GJ = " << TMath::Cos(locthetaphi[3]) << ", phi_omega_GJ = " << locthetaphi[4] << endl;
 
  cout << "Testing new utility class:" << endl;
  vector< double > locXTo4PiAnglesH = getTwoStepAngles(X, vec, vec_daught1, vec_daught2, pi0, beam, target, 1, true);
  vector< double > locXTo4PiAnglesGJ = getTwoStepAngles(X, vec, vec_daught1, vec_daught2, pi0, beam, target, 2, true);
  cout << "Helicity frame:" << endl;
  cout << "cos(theta_omega) = " << TMath::Cos(locXTo4PiAnglesH[0]) << ", phi_omega = " << locXTo4PiAnglesH[1] << ", cos(theta_H) = " << TMath::Cos(locXTo4PiAnglesH[2]) << ", phi_H = " << locXTo4PiAnglesH[3] << endl;
  cout << "Gottfried-Jackson frame:" << endl;
  cout << "cos(theta_omega) = " << TMath::Cos(locXTo4PiAnglesGJ[0]) << ", phi_omega = " << locXTo4PiAnglesGJ[1] << ", cos(theta_H) = " << TMath::Cos(locXTo4PiAnglesGJ[2]) << ", phi_H = " << locXTo4PiAnglesGJ[3] << endl;
  cout << "Lower Vertex Decay:" << endl;
  vector< double > locDeltaAnglesH = getOneStepAngles(deltaRecoil, recoil, beam, target, 1, false);
  vector< double > locDeltaAnglesGJ = getOneStepAngles(deltaRecoil, recoil, beam, target, 2, false);
  cout << "Helicity frame:" << endl;
  cout << "cos(theta_p) = " << TMath::Cos(locDeltaAnglesH[0]) << ", phi_p = " << locDeltaAnglesH[1] << endl;
  cout << "Gottfried-Jackson frame:" << endl;
  cout << "cos(theta_p) = " << TMath::Cos(locDeltaAnglesGJ[0]) << ", phi_p = " << locDeltaAnglesGJ[1] << endl;
*/
  return;
}


////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
Vec_ps_refl::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{

  GDouble cosTheta = userVars[uv_cosTheta];
  GDouble Phi = userVars[uv_Phi];
  GDouble cosThetaH = userVars[uv_cosThetaH];
  GDouble PhiH = userVars[uv_PhiH];
  GDouble prod_angle = userVars[uv_prod_Phi];
  GDouble MX = userVars[uv_MX];
  GDouble MVec = userVars[uv_MVec];
  GDouble MPs = userVars[uv_MPs];
  GDouble beam_polFraction = userVars[uv_beam_polFraction];
  GDouble beam_polAngle = userVars[uv_beam_polAngle];


  complex <GDouble> amplitude(0,0);
  complex <GDouble> i(0,1);

  for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
	  GDouble hel_amp = clebschGordan(m_l, 1, 0, lambda, m_j, lambda);
	  amplitude += conj(wignerD( m_j, m_m, lambda, cosTheta, Phi )) * hel_amp * conj(wignerD( 1, lambda, 0, cosThetaH, PhiH ));
  } 
  
  GDouble Factor = sqrt(1 + m_s * beam_polFraction);
  
  
  complex< GDouble > zjm = 0;
  
  complex< GDouble > rotateY = polar( (GDouble)1., (GDouble)(-1.*(prod_angle + beam_polAngle*TMath::DegToRad())) ); // - -> + in prod_angle and polAngle summing
  
  
  if (m_r == 1)
	  zjm = real(amplitude * rotateY);
  if (m_r == -1) 
	  zjm = i*imag(amplitude * rotateY);

  // E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
  double kinFactor = barrierFactor(MX, m_l, MVec, MPs);
  //kinFactor *= sqrt(3.) * sqrt(2.*m_l + 1.);
  Factor *= kinFactor;

  return complex< GDouble >( static_cast< GDouble>( Factor ) * zjm );
}


void Vec_ps_refl::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
Vec_ps_refl::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUVec_ps_refl_exec( dimGrid, dimBlock, GPU_AMP_ARGS, m_j, m_m, m_l, m_r, m_s );

}

#endif


