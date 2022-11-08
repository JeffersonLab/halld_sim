
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/omegapiAngles.h"

OmegaDalitz::OmegaDalitz( const vector< string >& args ) :
UserAmplitude< OmegaDalitz >( args )
{
	//assert( args.size() == 4 );
  
	//Dalitz Parameters for 3pi decays
	dalitz_alpha  = AmpParameter(args[0]);
	dalitz_beta   = AmpParameter(args[1]);
	dalitz_gamma  = AmpParameter(args[2]);
	dalitz_delta  = AmpParameter(args[3]);
	
	registerParameter(dalitz_alpha);
	registerParameter(dalitz_beta);
	registerParameter(dalitz_gamma);
	registerParameter(dalitz_delta);
}

void
OmegaDalitz::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 

  // common vector and pseudoscalar P4s
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); // 1st after proton
  TLorentzVector vec, vec_daught1, vec_daught2; // compute for each final state below 

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
  double dalitz_sc = (1/3.)*( vec.M2() - pip.M2() - pim.M2() - pi0.M2());
  double dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
  double dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;
  double dalitz_z = dalitzx*dalitzx + dalitzy*dalitzy;
  double dalitz_sin3theta = TMath::Sin(3 *  TMath::ASin( (dalitzy/sqrt(dalitz_z) )) );
  double dalitz_phi = dalitz_s*dalitz_t*dalitz_u - pi0.M2()*pow(vec.M2() - pi0.M2(), 2.);
  
  userVars[uv_dalitz_z] = dalitz_z;
  userVars[uv_dalitz_sin3theta] = dalitz_sin3theta;
  userVars[uv_dalitz_phi] = dalitz_phi;

  return;
}


////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
OmegaDalitz::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{

  GDouble dalitz_z = userVars[uv_dalitz_z];
  GDouble dalitz_sin3theta = userVars[uv_dalitz_sin3theta];
  GDouble dalitz_phi = userVars[uv_dalitz_phi];

  // dalitz parameters for 3-body vector decay
  GDouble G = sqrt( fabs(dalitz_phi * (1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * pow(dalitz_z,3/2.) * dalitz_sin3theta + 2 * dalitz_gamma * pow(dalitz_z,2) + 2 * dalitz_delta * pow(dalitz_z,5/2.) * dalitz_sin3theta)) );

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


