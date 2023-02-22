
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/omegapiAngles.h"
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

  polAngle = AmpParameter( args[5] );
  registerParameter( polAngle );
  
  polFraction = atof(args[6].c_str());
  
  // BeamProperties configuration file
  if (polFraction == 0){
    TString beamConfigFile = args[6].c_str();
    BeamProperties beamProp(beamConfigFile);
    polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
  }

  m_3pi = false;
  if(args.size() == (11)){
	  m_3pi = true; // treat 3-pion decay dalitz parameters

	  //Dalitz Parameters for 3pi decays
	  dalitz_alpha  = AmpParameter(args[6+1]);
	  dalitz_beta   = AmpParameter(args[6+2]);
	  dalitz_gamma  = AmpParameter(args[6+3]);
	  dalitz_delta  = AmpParameter(args[6+4]);
	  
	  registerParameter(dalitz_alpha);
	  registerParameter(dalitz_beta);
	  registerParameter(dalitz_gamma);
	  registerParameter(dalitz_delta);
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
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
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
	  
	  ///////////////////////////////////////////// Dalitz Parameters ///////////////////////////////
	  double dalitz_s = (pip+pim).M2(); //s=M2(pip pim)
	  double dalitz_t = (pip+pi0).M2(); //t=M2(pip pi0)
	  double dalitz_u = (pim+pi0).M2(); //u=M2(pim pi0)
	  double m3pi = (2*pip.M())+pi0.M();
	  double dalitz_d = 2*vec.M()*( vec.M() - m3pi);
	  double dalitz_sc = (1/3.)*( vec.M2() + pip.M2() + pim.M2() + pi0.M2());
	  double dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
	  double dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;
	  double dalitz_z = dalitzx*dalitzx + dalitzy*dalitzy;
	  double dalitz_sin3theta = TMath::Sin(3 *  TMath::ASin( (dalitzy/sqrt(dalitz_z) )) );

	  userVars[uv_dalitz_z] = dalitz_z;
	  userVars[uv_dalitz_sin3theta] = dalitz_sin3theta;
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
  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;

  // Calculate decay angles in helicity frame (same for all vectors)
  // set beam polarization angle to 0 degrees; apply diamond orientation in calcAmplitude
  vector <double> locthetaphi = getomegapiAngles(0, vec, X, beam, Gammap);

  // Calculate vector decay angles (unique for each vector)
  vector <double> locthetaphih;
  if(m_3pi) locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, vec_daught2);
  else locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, TLorentzVector(0,0,0,0));

  userVars[uv_cosTheta] = TMath::Cos(locthetaphi[0]);
  userVars[uv_Phi] = locthetaphi[1];

  userVars[uv_cosThetaH] = TMath::Cos(locthetaphih[0]);
  userVars[uv_PhiH] = locthetaphih[1];

  userVars[uv_prod_Phi] = locthetaphi[2];
  userVars[uv_dalitz_phi] = locthetaphih[2];

  userVars[uv_MX] = X.M();
  userVars[uv_MVec] = vec.M();
  userVars[uv_MPs] = ps.M();

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
  GDouble dalitz_z = userVars[uv_dalitz_z];
  GDouble dalitz_sin3theta = userVars[uv_dalitz_sin3theta];
  GDouble dalitz_phi = userVars[uv_dalitz_phi];
  GDouble MX = userVars[uv_MX];
  GDouble MVec = userVars[uv_MVec];
  GDouble MPs = userVars[uv_MPs];

  // dalitz parameters for 3-body vector decay
  GDouble G = 1; // not relevant for 2-body vector decays 
  if(m_3pi) G = sqrt( fabs(dalitz_phi * (1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * pow(dalitz_z,3/2.) * dalitz_sin3theta + 2 * dalitz_gamma * pow(dalitz_z,2) + 2 * dalitz_delta * pow(dalitz_z,5/2.) * dalitz_sin3theta)) );

  complex <GDouble> amplitude(0,0);
  complex <GDouble> i(0,1);

  for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
	  GDouble hel_amp = clebschGordan(m_l, 1, 0, lambda, m_j, lambda);
	  amplitude += conj(wignerD( m_j, m_m, lambda, cosTheta, Phi )) * hel_amp * conj(wignerD( 1, lambda, 0, cosThetaH, PhiH )) * G;
  } 

  GDouble Factor = sqrt(1 + m_s * polFraction);
  complex< GDouble > zjm = 0;
  complex< GDouble > rotateY = polar( (GDouble)1., (GDouble)(-1.*(prod_angle + polAngle*TMath::DegToRad())) ); // - -> + in prod_angle and polAngle summing
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

	GPUVec_ps_refl_exec( dimGrid, dimBlock, GPU_AMP_ARGS, m_j, m_m, m_l, m_r, m_s, m_3pi, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction );

}

#endif


