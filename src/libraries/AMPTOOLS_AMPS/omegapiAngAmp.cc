//July 26th 2019, Based on DOI: 10.1016/0550-3213(84)90382-1
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include "UTILITIES/CobremsGeneration.hh"
#include "UTILITIES/BeamProperties.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "AMPTOOLS_AMPS/omegapiAngles.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

TLorentzVector targetopi(0,0,0,0.938);

 double pararray[22];

//Create array of lmLM:
int lmLM[25][4] = {{0,0,0,0}, {0,0,2,0}, {0,0,2,1}, {0,0,2,2}, {2,0,0,0}, {2,0,2,0}, {2,0,2,1}, {2,0,2,2}, {2,1,2,0}, {2,1,2,1}, {2,1,2,2}, {2,2,2,0}, {2,2,2,1}, {2,2,2,2}, {2,1,1,1}, {0,0,1,0}, {0,0,1,1}, {2,1,1,0}, {2,1,1,1}, {2,1,2,1}, {2,1,2,2}, {2,2,2,1}, {2,2,2,2}, {2,0,1,0}, {2,0,1,1}};

 
int delta(int first, int second)
{
  int kroneckerdelta = 0;

  if (first == second) kroneckerdelta = 1;
    
  //cout << "kroneker delta =" << kroneckerdelta << endl;
  return kroneckerdelta;
}

double calpha(int alpha)
{
  double normalization = 0.0;
  
  int l = lmLM[alpha][0];
  int m = lmLM[alpha][1];
  int L = lmLM[alpha][2];
  int M = lmLM[alpha][3];
  
  double normalization_num = 16.0 * TMath::Pi() * TMath::Pi();
  double normalization_denum = ((2.0 * l)+1.0) * ((2.0 * L)+1.0) * (2.0-delta(m,0)) * (2.0-delta(M,0));
  normalization = normalization_denum/normalization_num;
  
  return normalization;
}

double hmoment(int alpha, vector<double> vector)
{
  double loccostheta = TMath::Cos(vector[0]);
  double locphi = vector[1];
  double loccosthetaH = TMath::Cos(vector[2]);
  double locphiH = vector[3];

  int l = lmLM[alpha][0];
  int m = lmLM[alpha][1];
  int L = lmLM[alpha][2];
  int M = lmLM[alpha][3];

   double moment = 0.0;
  if (alpha < 15)
  {
    moment = 0.5 * std::real(wignerD(L, M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH) + pow(-1.0,L+M) * wignerD(L, -M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH));
  }
  else {moment = 0.5 * std::real(wignerD(L, M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH) - pow(-1.0,L+M) * wignerD(L, -M, m, loccostheta, locphi) * wignerD(l, m, 0, loccosthetaH, locphiH));}

  return moment;
}

omegapiAngAmp::omegapiAngAmp( const vector< string >& args ):
  UserAmplitude< omegapiAngAmp >( args )
{
	assert( args.size() == 24 || args.size() == 25 );
	
	if(args.size() == 25){
		polAngle  = atof(args[22].c_str() ); // azimuthal angle of the photon polarization vector in the lab measured in degrees.
		polFraction = AmpParameter( args[23] ); // polarization fraction
		std::cout << "Fixed polarization fraction =" << polFraction << " and pol.angle= " << polAngle << " degrees." << std::endl;
	}
	else if (args.size() == 24){
		// BeamProperties configuration file
		TString beamConfigFile = args[22].c_str();
		BeamProperties beamProp(beamConfigFile);
		polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
		polAngle = beamProp.GetPolAngle();
		std::cout << "Polarisation angle of " << polAngle << " from BeamProperties." << std::endl;
		if(polAngle == -1)
			std::cout << "This is an amorphous run. Set beam polarisation to 0." << std::endl;
		for(Int_t i=0; i<polFrac_vs_E->GetXaxis()->GetNbins()+2; i++){
			//cout << polFrac_vs_E->GetBinContent(i) << endl;
		}
	}
	else
	assert(0);

    //A switch to use cut off only when generating (casues problem when fitting)
    //Set to 1 with gen_amp, 0 with fit.
    useCutoff = atoi( args[24].c_str() );
 
	
    ds_ratio = AmpParameter(args[9]);
    
    //1+ state
    m_1p = AmpParameter(args[0]);
    w_1p = AmpParameter(args[1]);
    n_1p = AmpParameter(args[2]);
    phi0_1p = AmpParameter(args[10]);
    theta_1p = AmpParameter(args[11]);
    phip_1p = AmpParameter(args[12]);
    phim_1p = AmpParameter(args[13]);
    psi_1p = AmpParameter(args[14]);
    //1- state    
    m_1m = AmpParameter(args[3]);
    w_1m = AmpParameter(args[4]);
    n_1m = AmpParameter(args[5]);
    phi0_1m = AmpParameter(args[15]);
    theta_1m = AmpParameter(args[16]);
    phip_1m = AmpParameter(args[17]);
    phim_1m = AmpParameter(args[18]);
    psi_1m = AmpParameter(args[19]);
    //0- state    
    m_0m = AmpParameter(args[6]);
    w_0m = AmpParameter(args[7]);
    n_0m = AmpParameter(args[8]);
    phi0_0m = AmpParameter(args[20]);
    theta_0m = AmpParameter(args[21]);    

  registerParameter(m_1p);
  registerParameter(w_1p);
  registerParameter(n_1p);
  registerParameter(phi0_1p);
  registerParameter(theta_1p);
  registerParameter(phip_1p);
  registerParameter(phim_1p);
  registerParameter(psi_1p);

  registerParameter(m_1m);
  registerParameter(w_1m);
  registerParameter(n_1m);
  registerParameter(phi0_1m);
  registerParameter(theta_1m);
  registerParameter(phip_1m);
  registerParameter(phim_1p);
  registerParameter(psi_1m);

  registerParameter(m_0m);
  registerParameter(w_0m);
  registerParameter(n_0m);
  registerParameter(phi0_0m);
  registerParameter(theta_0m);

  registerParameter(ds_ratio);

    pararray[0] = m_1p;
    pararray[1] = w_1p;
    pararray[2] = n_1p;
    pararray[3] = m_1m;
    pararray[4] = w_1m;
    pararray[5] = n_1m;
    pararray[6] = m_0m;
    pararray[7] = w_0m;
    pararray[8] = n_0m;
    pararray[9] = ds_ratio;
    pararray[10] = phi0_1p;
    pararray[11] = theta_1p;
    pararray[12] = phip_1p;
    pararray[13] = phim_1p;
    pararray[14] = psi_1p;
    pararray[15] = phi0_1m;
    pararray[16] = theta_1m;
    pararray[17] = phip_1m;
    pararray[18] = phim_1m;
    pararray[19] = psi_1m;
    pararray[20] = phi0_0m;
    pararray[21] = theta_0m;
}

//Define breakup momentum (x here is the measured mass of the omegapi)
double q1(double x) {
  double pi0mass = 0.1349766; //mass of the pi0 in GeV
  double omegamass = 0.78265; //mass of the omega in GeV
  if (x < pi0mass + omegamass)
    return 0.0;
  return 0.5 * TMath::Sqrt((TMath::Power(x, 2.0) - (pi0mass + omegamass)*(pi0mass + omegamass))*(TMath::Power(x, 2.0) - (omegamass - pi0mass)*(omegamass - pi0mass))) / x;
}

//Define barrier factor ratio
double barrierratio(double x, double resonancemass, int l) {
  return barrierFactor(q1(x), l) / barrierFactor(q1(resonancemass), l);
}

double Gamma_alpha(double x, double resonancewidth, double resonancemass, int l) {
  return resonancewidth * q1(x) / q1(resonancemass) * TMath::Power(barrierratio(x, resonancemass, l), 2);
}

std::complex<double> D_alpha(double x, double resonancemass, double resonancewidth, int l) {
  std::complex<double> denom = std::complex<double>(TMath::Power(x, 2) - TMath::Power(resonancemass, 2),-1.0*resonancemass * Gamma_alpha(x, resonancewidth, resonancemass, l));
  return resonancewidth * resonancemass / denom;
}

//Define J (spin) for a given ialpha
int J_spin(int ialpha) { //ialpha{0,1,2} -> JP{1+,1-,0-}
  if (ialpha == 2)
    return 0;
  return 1;
}

//Define a parity function
int eta_parity(int ialpha) {
  if (ialpha == 0)
    return 1;
  return -1;
}

//Define sum over l with Clebsch-Gordan coefficients and barrier factors
double lsum(double x, double resonancemass, int ialpha, int lambda, double DoverS) {
  double c_alpha; //partial wave amplitudes
  double lsum = 0;
  for (int l = 0; l < 3; l++) { 
    if (l == 1 && (ialpha == 1 || ialpha == 2))
      c_alpha = 1;
    else if (l == 2 && ialpha == 0)
      c_alpha = DoverS;
    else if (l == 0 && ialpha == 0)
      c_alpha = 1;
    else 
      c_alpha = 0; 
    lsum += TMath::Sqrt((2.0*l + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * clebschGordan(l, 1, 0, lambda, J_spin(ialpha), lambda) * c_alpha * barrierratio(x, resonancemass, l);
  }
  return lsum;
}

std::complex<double> F_alpha_lambda(double x, double resonancemass, double resonancewidth, double G_alpha, int ialpha, int lambda, double DoverS, int l) {
  return D_alpha(x, resonancemass, resonancewidth, l) * G_alpha * lsum(x, resonancemass, ialpha, lambda, DoverS);
}

//Define sum over omega spin states
std::complex<double> f_Llm(double x, double resonancemass, double resonancewidth, double G_alpha, int ialpha, double resonancemass2, double resonancewidth2, double G_beta, int ibeta, double DoverS, int l, int m, int L) {
  std::complex<double> fsum = std::complex<double>(0, 0);
  for (int lambda = -1; lambda < 2; lambda++) {
    for (int lambdaprime = -1; lambdaprime < 2; lambdaprime++) {
      if (ibeta == 2 && lambdaprime != 0)
	continue;
      if (ialpha == 2 && lambda != 0)
	continue;
      fsum += F_alpha_lambda(x, resonancemass, resonancewidth, G_alpha, ialpha, lambda, DoverS, l) * std::conj(F_alpha_lambda(x, resonancemass2, resonancewidth2, G_beta, ibeta, lambdaprime, DoverS, l)) * clebschGordan(J_spin(ibeta), L, lambdaprime, m, J_spin(ialpha), lambda) * clebschGordan(1, l, lambdaprime, m, 1, lambda);
    }
  }
  return fsum;
}

const std::complex<double> ic(0, 1);

//Define complex "amplitudes"
std::complex<double> f_0(double phi0, double theta) { 
  return TMath::Sqrt(0.5) * std::exp(ic * phi0) * TMath::Cos(theta);
}
std::complex<double> f_plus(double phiplus, double theta, double psi) {
  return TMath::Sqrt(0.5) * std::exp(ic * phiplus) * TMath::Sin(theta) * TMath::Cos(psi);
}
std::complex<double> f_minus(double phiminus, double theta, double psi) {
  return TMath::Sqrt(0.5) * std::exp(ic * phiminus) * TMath::Sin(theta) * TMath::Sin(psi);
}

//Define a helicity function
int Lambda_H(int iH) { // {0,1,2}->{0,+1,-1}
  if (iH == 0)
    return 0;
  else if (iH == 1) 
    return 1;
  else
    return -1;
}


//Define production density matrix
std::complex<double> rho(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int iH, int iHPrime, int ialpha, int ibeta) {
  std::complex<double> f_alpha = std::complex<double>(0,0);
  std::complex<double> f_beta = std::complex<double>(0,0); 
  std::complex<double> f_alpha_neg = std::complex<double>(0,0);
  std::complex<double> f_beta_neg = std::complex<double>(0,0); 
  if (iH == 0) {
    f_alpha = f_0(phi0, theta);
    f_alpha_neg = f_0(phi0, theta);
  }
  else if (iH == 1) {
    f_alpha = f_plus(phiplus, theta, psi);
    f_alpha_neg = f_minus(phiminus, theta, psi);
  }
  else { 
    f_alpha = f_minus(phiminus, theta, psi);
    f_alpha_neg = f_plus(phiplus, theta, psi);
  }
  if (iHPrime == 0) {
    f_beta = f_0(phi02, theta2);
    f_beta_neg = f_0(phi02, theta2);
  }
  else if (iHPrime == 1) {
    f_beta = f_plus(phiplus2, theta2, psi2);
    f_beta_neg = f_minus(phiminus2, theta2, psi2);
  }
  else { 
    f_beta = f_minus(phiminus2, theta2, psi2);
    f_beta_neg = f_plus(phiplus2, theta2, psi2);
  }

  return f_alpha * std::conj(f_beta) + eta_parity(ialpha) * eta_parity(ibeta) * TMath::Power(-1.0, J_spin(ialpha) - J_spin(ibeta)) * TMath::Power(-1.0, Lambda_H(iH) - Lambda_H(iHPrime)) * f_alpha_neg * std::conj(f_beta_neg);

}


//Define sum over helicity states
std::complex<double> HelicitySum(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int ialpha, int ibeta, int L, int M) {
  std::complex<double> sumH(0, 0);
  for (int iH = 0; iH < 3; iH++) {
    for (int iHPrime = 0; iHPrime < 3; iHPrime++) {
      if (ialpha == 2 && iH > 0)
	continue;
      if (ibeta == 2 && iHPrime > 0)
	continue;
      sumH += rho(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, iH, iHPrime, ialpha, ibeta) * clebschGordan(J_spin(ibeta), L, Lambda_H(iHPrime), M, J_spin(ialpha), Lambda_H(iH));
    }
  }
  
      //cout << "SumH = " << sumH << endl;
  return sumH;
}

//Define t_star (sum of production density matrix)
std::complex<double> t_star_LM(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int ialpha, int ibeta, int L, int M) {

  complex <GDouble> t_star_LM_par = TMath::Sqrt((2.0*J_spin(ibeta) + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * HelicitySum(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, ialpha, ibeta, L, M);
  
      //cout << "T*(LM) = " << t_star_LM_par << endl;

  return t_star_LM_par;
}

double SingleIntensity(double x, double *pararray, int l, int m, int L, int M, int ialpha, int ibeta) {
    double single_intensity = 0.;
    double phiplus_alpha = 0;
    double phiminus_alpha = 0;
    double psi_alpha = 0;
    double phiplus_beta = 0;
    double phiminus_beta = 0;
    double psi_beta = 0;
    if (ialpha < 2) {
      phiplus_alpha = pararray[5*ialpha + 12];
      phiminus_alpha = pararray[5*ialpha + 13];
      psi_alpha = pararray[5*ialpha + 14];
    }
    if (ibeta < 2) {
      phiplus_beta = pararray[5*ibeta + 12];
      phiminus_beta = pararray[5*ibeta + 13];
      psi_beta = pararray[5*ibeta + 14];
    }

    single_intensity = std::real(t_star_LM(pararray[5*ialpha + 10], pararray[5*ialpha + 11], phiplus_alpha, phiminus_alpha, psi_alpha, pararray[5*ibeta + 10], pararray[5*ibeta + 11], phiplus_beta, phiminus_beta, psi_beta, ialpha, ibeta, L, M) * f_Llm(x, pararray[3*ialpha + 0], pararray[3*ialpha + 1], pararray[3*ialpha + 2], ialpha, pararray[3*ibeta + 0], pararray[3*ibeta + 1], pararray[3*ibeta + 2], ibeta, pararray[9], l, m, L) * clebschGordan(1, l, 0, 0, 1, 0));
    //cout << "single intensity = " << single_intensity << endl;

  return single_intensity;
}

double Intensity(double x, double *pararray, int i){
  
    double IntensitySum = 0;
  for (int ialpha = 0; ialpha < 3; ialpha++) { //double sum over combination of JP states
    for (int ibeta = 0; ibeta < 3; ibeta++) {
      if ( (ialpha != ibeta) & (eta_parity(ialpha) == eta_parity(ibeta)) ) //Only want interference moments with opposite parity
	continue;
      IntensitySum += SingleIntensity(x, pararray, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
    }
  }
    //cout << "IntensitySum = " << IntensitySum << endl;

  return IntensitySum;
}

double wdist(int k, double x, double Phi, vector<double> vector)
  {
    double dist = 0.0;

          for (int alpha = 0; alpha < 25; alpha++)//quantum states
	  {
	    if (k == 0){dist += Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha);}

	    if (k == 1){dist += Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha) * TMath::Sqrt(2.0) * TMath::Cos(2.0 * Phi);}

	    if (k == 2){dist += Intensity(x, pararray, alpha) * hmoment(alpha, vector) * calpha(alpha) * TMath::Sqrt(2.0) * TMath::Sin(2.0 * Phi);}

	  }//alpha loop
  // cout << "wdist(" << k << ") = " << dist << endl;

    return dist;
  }

    complex <GDouble> sqrtIntensity(double polfrac, double x, double Phi, vector<double> vector)
//    double sqrtIntensity(double polfrac, double x, double Phi, vector<double> vector)
  {
	double intensity = 0.0;
	double real_sqrt_intensity = 0.0;
	double ima_sqrt_intensity = 0.0;

	intensity = wdist(0, x, Phi,vector) - polfrac * wdist(1, x, Phi,vector) * TMath::Cos(2.0 * Phi) -  polfrac * wdist(2, x, Phi,vector) * TMath::Sin(2.0 * Phi);
	
       if (intensity >= 0.0) real_sqrt_intensity = sqrt(intensity);
       if (intensity < 0.0) ima_sqrt_intensity = sqrt(abs(intensity));
	
	complex <GDouble> sqrt_intensity(real_sqrt_intensity,ima_sqrt_intensity);
	//cout << "intensity = " << intensity << ", sqrt_intensity = " << sqrt_intensity << endl;
       return sqrt_intensity;
  }

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
omegapiAngAmp::calcAmplitude( GDouble** pKin ) const
{
  complex <GDouble> CZero(0,0);  

  ////////////////////////////////////////////////// Get Lorentz Vectors of Particles //////////////////////////////////

  TLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  TLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

  TLorentzVector rhos_pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector rhos_pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
  TLorentzVector rho = rhos_pip + rhos_pim;

  TLorentzVector omegas_pi(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector omega = rho + omegas_pi;

  //Cut off to force a narrow omega peak. Casues problem with fit.
  GDouble m0_omega=0.782;
  GDouble mG0_omega = 0.00849;    
  if(useCutoff && fabs(omega.M()-m0_omega) > mG0_omega){
    //cout << "Cutting event with m_omega =" << omega.M() << endl;
    return CZero; 
  }

  TLorentzVector Xs_pi(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);
  TLorentzVector X = omega + Xs_pi;

    ////////////////////////////////////////////////// Boost Particles and Get Angles//////////////////////////////////

  //Helicity coordinate system
  TLorentzVector Gammap = beam + targetopi;
 
// polarization BeamProperties
	GDouble Pgamma=polFraction;//fixed beam polarization fraction
	if(polAngle == -1)
	Pgamma = 0.;//if beam is amorphous set polarization fraction to 0
	else if(polFrac_vs_E!=NULL){
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());

	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else
	 Pgamma = polFrac_vs_E->GetBinContent(bin);
	}
   double mx = X.M();

  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);
  
  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

  //cout << "theta =" << locthetaphi[0] << ", phi= " << locthetaphi[1] << ", thetah =" << locthetaphih[0] << ", phih= " << locthetaphih[1] << endl;
  vector <double> angvector{locthetaphi[0], locthetaphi[1], locthetaphih[0], locthetaphih[1]};

    GDouble Phi = locthetaphi[2];

    complex <GDouble> amplitude = sqrtIntensity(Pgamma, mx, Phi, angvector);
 
  //cout << real(amplitude) << "+i" << imag(amplitude) << endl;


return amplitude;
}

void omegapiAngAmp::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates  
}
