//Jan 18th 2020, Based on model by Adam Szczepaniak & Vincent Mathieu
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
//#include "UTILITIES/CobremsGeneration.hh"
//#include "UTILITIES/BeamProperties.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "omegapi_amplitude.h"
#include "barrierFactor.h"
#include "clebschGordan.h"
#include "wignerD.h"
#include "breakupMomentum.h"
#include "omegapiAngles.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"


omegapi_amplitude::omegapi_amplitude( const vector< string >& args ):
  UserAmplitude< omegapi_amplitude >( args )
{
	assert( args.size() == (7+4+2) || args.size() == (7+4+3) );
	
	if(args.size() == (7+4+3)){
		polAngle  = atof(args[7+4+1].c_str() ); // azimuthal angle of the photon polarization vector in the lab measured in degrees.
		polFraction = AmpParameter( args[7+4+2] ); // polarization fraction
		std::cout << "Fixed polarization fraction =" << polFraction << " and pol.angle= " << polAngle << " degrees." << std::endl;
	}
/*	else if (args.size() == (7+4+2)){//beam properties requires halld_sim
		// BeamProperties configuration file
		TString beamConfigFile = args[7+4+1].c_str();
		BeamProperties beamProp(beamConfigFile);
		polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
		polAngle = beamProp.GetPolAngle();
		std::cout << "Polarisation angle of " << polAngle << " from BeamProperties." << std::endl;
		if(polAngle == -1)
			std::cout << "This is an amorphous run. Set beam polarisation to 0." << std::endl;
		for(Int_t i=0; i<polFrac_vs_E->GetXaxis()->GetNbins()+2; i++){
			//cout << polFrac_vs_E->GetBinContent(i) << endl;
		}
	}*/
	else
	assert(0);

    sign = atoi(args[0].c_str() );
    lambda_gamma = atoi(args[1].c_str() );
    spin = atoi(args[2].c_str() );
    parity = atoi(args[3].c_str() );
    spin_proj = atoi(args[4].c_str() );
    
    c_0 = AmpParameter(args[5]);
    registerParameter(c_0);
    
    c_1 = AmpParameter(args[6]);
    registerParameter(c_1);
    
    c_2 = AmpParameter(args[7]);
    registerParameter(c_2);
    
   //Dalitz Parameters
   dalitz_alpha  = AmpParameter(args[7+1]);
   dalitz_beta  = AmpParameter(args[7+2]);
   dalitz_gamma  = AmpParameter(args[7+3]);
   dalitz_delta  = AmpParameter(args[7+4]);

   registerParameter(dalitz_alpha);
   registerParameter(dalitz_beta);
   registerParameter(dalitz_gamma);
   registerParameter(dalitz_delta);

}
////////////////////////////////////////////////// User Vars //////////////////////////////////
void
omegapi_amplitude::calcUserVars( GDouble** pKin, GDouble* userVars ) const 
{

  TLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  TLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

  TLorentzVector rhos_pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector rhos_pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
  TLorentzVector rho = rhos_pip + rhos_pim;

  TLorentzVector omegas_pi(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector omega = rho + omegas_pi;

  TLorentzVector Xs_pi(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);

  TLorentzVector X = omega + Xs_pi;

    //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  TLorentzVector target(0,0,0,0.938);
  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;
 
// polarization BeamProperties
	GDouble Pgamma=polFraction;//fixed beam polarization fraction
	if(polAngle == -1)
	Pgamma = 0.;//if beam is amorphous set polarization fraction to 0
/*	else if(polFrac_vs_E!=NULL){
	//This part causes seg fault with 34 amplitudes or more with gen_amp and gen_omegapi.
	//Not needed for fixed beam pol angle and frac.
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());

	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else
	 Pgamma = polFrac_vs_E->GetBinContent(bin);
	}*/

  //Calculate decay angles in helicity frame
  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

  userVars[uv_cosTheta] = TMath::Cos(locthetaphi[0]);
  userVars[uv_Phi] = locthetaphi[1];

  userVars[uv_cosThetaH] = TMath::Cos(locthetaphih[0]);
  userVars[uv_PhiH] = locthetaphih[1];

  userVars[uv_prod_angle] = locthetaphi[2];

  userVars[uv_Pgamma] = Pgamma;
  
///////////////////////////////////////////// Dalitz Parameters ///////////////////////////////
  double dalitz_s = rho.M2();//s=M2(pip pim)
  double dalitz_t = (rhos_pip+omegas_pi).M2();//t=M2(pip pi0)
  double dalitz_u = (rhos_pim+omegas_pi).M2();//u=M2(pim pi0)
  double m3pi = (2*139.57018)+134.9766;
  double dalitz_d = 2*omega.M()*( omega.M() - m3pi);
  double dalitz_sc = (1/3)*( omega.M2() - rhos_pip.M2() - rhos_pim.M2() - omegas_pi.M2());
  double dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
  double dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;
  double dalitz_z = dalitzx*dalitzx + dalitzy*dalitzy;
  double dalitz_sin3theta = TMath::Sin(3 *  TMath::ASin( (dalitzy/sqrt(dalitz_z) )) );
  
  userVars[uv_dalitz_z] = dalitz_z;
  userVars[uv_dalitz_sin3theta] = dalitz_sin3theta;
  
}

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
omegapi_amplitude::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
    
  complex <GDouble> COne(1,0);  
  complex <GDouble> CZero(0,0);  

   GDouble cosTheta = userVars[uv_cosTheta];
   GDouble Phi = userVars[uv_Phi];
   GDouble cosThetaH = userVars[uv_cosThetaH];
   GDouble PhiH = userVars[uv_PhiH];
   GDouble prod_angle = userVars[uv_prod_angle];
   GDouble polfrac = userVars[uv_Pgamma];
   GDouble dalitz_z = userVars[uv_dalitz_z];
   GDouble dalitz_sin3theta = userVars[uv_dalitz_sin3theta];

   GDouble G = sqrt(1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * pow(dalitz_z,3/2) * dalitz_sin3theta
			 + 2 * dalitz_gamma * pow(dalitz_z,2) + 2 * dalitz_delta * pow(dalitz_z,5/2) * dalitz_sin3theta );

   GDouble hel_c[3] = { c_0, c_1, c_2};
   
 complex <GDouble> amplitude = CZero;
 
   for (int lambda = -1; lambda <= 1; lambda++)//omega helicity
	      {
		  GDouble hel_amp = 0.0;

		  for(int l = 0; l <= 2; l++)//partial waves (l).
		  {//if ( (parity == -1 && l% 2 == 0) || (parity == 1 && l%2 != 0) ) continue;
		  
		  hel_amp += hel_c[l] * clebschGordan(l, 1, 0, lambda, spin, lambda);
		  }//loop over l
		  
		  amplitude += wignerD( spin, spin_proj, lambda, cosTheta, Phi ) * hel_amp * wignerD( 1, lambda, 0, cosThetaH, PhiH ) * G;
		}//loop over lambda
		

 complex <GDouble> prefactor ( cos( lambda_gamma * prod_angle ), sin( lambda_gamma * prod_angle ));
 
 if (sign == -1 && lambda_gamma == -1){amplitude *= -1*sqrt( ( 1 + (sign * polfrac) )/2 ) * prefactor;}
 else{amplitude *= sqrt( ( 1 + (sign * polfrac) )/2 ) * prefactor;}

//Instead of the vertices "scale" in configuration file.
// if( (parity == +1 && ((spin+spin_proj)%2 == 0) ) || (parity == -1 && ((spin+spin_proj)%2 != 0) ) ){amplitude *= -1;} 

return amplitude;
}

void omegapi_amplitude::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates  
}

#ifdef GPU_ACCELERATION
void
omegapi_amplitude::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
    
  GPUomegapi_amplitude_exec( dimGrid, dimBlock, GPU_AMP_ARGS, 
			  sign, lambda_gamma, spin, parity, spin_proj, c_0, c_1, c_2, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction);
}
#endif //GPU_ACCELERATION

