
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/vecPsAngles.h"

OmegaDalitz::OmegaDalitz( const vector< string >& args ) :
UserAmplitude< OmegaDalitz >( args )
{
	//assert( args.size() == 4 );
  
	// Dalitz Parameters for 3 pion decays
	dalitz_alpha  = AmpParameter(args[0]);
	dalitz_beta   = AmpParameter(args[1]);
	dalitz_gamma  = AmpParameter(args[2]);
	dalitz_delta  = AmpParameter(args[3]);
	
	registerParameter(dalitz_alpha);
	registerParameter(dalitz_beta);
	registerParameter(dalitz_gamma);
	registerParameter(dalitz_delta);

	// index in event kinematics array for 3 pions (default to match Vec_ps_refl after bachelor PS) 
	index_daughter1 = 3;
	index_daughter2 = 4;
	index_daughter3 = 5;
	if(args.size() > 4) { // user provides index of omega decay daughters
		index_daughter1 = atoi(args[4].c_str());
		index_daughter2 = atoi(args[5].c_str());
		index_daughter3 = atoi(args[6].c_str());
	}
}

void
OmegaDalitz::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  // 
  TLorentzVector pi0(pKin[index_daughter1][1], pKin[index_daughter1][2], pKin[index_daughter1][3], pKin[index_daughter1][0]);
  TLorentzVector pip(pKin[index_daughter2][1], pKin[index_daughter2][2], pKin[index_daughter2][3], pKin[index_daughter2][0]);
  TLorentzVector pim(pKin[index_daughter3][1], pKin[index_daughter3][2], pKin[index_daughter3][3], pKin[index_daughter3][0]);
  TLorentzVector omega = pi0 + pip + pim;
  
  ///////////////////////////////////////////// Dalitz Parameters ///////////////////////////////
  double dalitz_s = (pip+pim).M2(); //s=M2(pip pim)
  double dalitz_t = (pip+pi0).M2(); //t=M2(pip pi0)
  double dalitz_u = (pim+pi0).M2(); //u=M2(pim pi0)
  double m3pi = (2*pip.M())+pi0.M();
  double dalitz_d = 2*omega.M()*( omega.M() - m3pi);
  double dalitz_sc = (1/3.)*( omega.M2() + pip.M2() + pim.M2() + pi0.M2());
  double dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
  double dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;
  double dalitz_z = dalitzx*dalitzx + dalitzy*dalitzy;
  double dalitz_sin3theta = TMath::Sin(3 *  TMath::ASin( (dalitzy/sqrt(dalitz_z) )) );

  TVector3 omegaboost = omega.BoostVector();
  pip.Boost(-1.0*omegaboost);
  pim.Boost(-1.0*omegaboost);
  TVector3 pip_omega = pip.Vect();
  TVector3 pim_omega = pim.Vect();
  TVector3 piCross = (pip_omega).Cross(pim_omega);
  double lambda = 4/3. * fabs(piCross.Dot(piCross)) / TMath::Power(1/9. * (omega.M2() - TMath::Power(2*pip.M() + pi0.M(), 2.)), 2.);
  
  userVars[uv_lambda] = lambda;

  userVars[uv_alpha_term] = 2 * dalitz_z;
  userVars[uv_beta_term] = 2 * pow(dalitz_z,3/2.) * dalitz_sin3theta;
  userVars[uv_gamma_term] = 2 * pow(dalitz_z,2);
  userVars[uv_delta_term] = 2 * pow(dalitz_z,5/2.) * dalitz_sin3theta;

  return;
}


////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
OmegaDalitz::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{

  GDouble lambda = userVars[uv_lambda];
  GDouble alpha_term = userVars[uv_alpha_term];
  GDouble beta_term = userVars[uv_beta_term];
  GDouble gamma_term = userVars[uv_gamma_term];
  GDouble delta_term = userVars[uv_delta_term];
  
  // dalitz parameters for 3-body omega decay
  
  GDouble G = sqrt( fabs( lambda *  ( 1 + dalitz_alpha * alpha_term + dalitz_beta * beta_term + dalitz_gamma * gamma_term + dalitz_delta * delta_term ) ) );
  
  return complex< GDouble >( G );
}

void OmegaDalitz::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
OmegaDalitz::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUOmegaDalitz_exec( dimGrid, dimBlock, GPU_AMP_ARGS, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta );

}

#endif


