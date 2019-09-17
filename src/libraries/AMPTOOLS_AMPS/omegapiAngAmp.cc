//July 26th 2019, Based on DOI: 10.1016/0550-3213(84)90382-1
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
// #include "UTILITIES/CobremsGeneration.hh"
// #include "UTILITIES/BeamProperties.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "omegapiAngAmp.h"
#include "barrierFactor.h"
#include "clebschGordan.h"
#include "wignerD.h"
#include "breakupMomentum.h"
#include "omegapiAngles.h"

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
		polAngle  = atof(args[23].c_str() ); // azimuthal angle of the photon polarization vector in the lab measured in degrees.
		polFraction = AmpParameter( args[24] ); // polarization fraction
		std::cout << "Fixed polarization fraction =" << polFraction << " and pol.angle= " << polAngle << " degrees." << std::endl;
	}
	else if (args.size() == 24){
		// BeamProperties configuration file
		TString beamConfigFile = args[23].c_str();
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
    //Set to 1 with gen_amp, 0 with fit. Not needed with the two-step generator.
    useCutoff = atoi( args[22].c_str() );
 
	
    m_ds_ratio = AmpParameter(args[9]);
    
    //1+ state
    m_1p = AmpParameter(args[0]);
    m_w_1p = AmpParameter(args[1]);
    m_n_1p = AmpParameter(args[2]);
    m_phi0_1p = AmpParameter(args[10]);
    m_theta_1p = AmpParameter(args[11]);
    m_phip_1p = AmpParameter(args[12]);
    m_phim_1p = AmpParameter(args[13]);
    m_psi_1p = AmpParameter(args[14]);
    //1- state    
    m_1m = AmpParameter(args[3]);
    m_w_1m = AmpParameter(args[4]);
    m_n_1m = AmpParameter(args[5]);
    m_phi0_1m = AmpParameter(args[15]);
    m_theta_1m = AmpParameter(args[16]);
    m_phip_1m = AmpParameter(args[17]);
    m_phim_1m = AmpParameter(args[18]);
    m_psi_1m = AmpParameter(args[19]);
    //0- state    
    m_0m = AmpParameter(args[6]);
    m_w_0m = AmpParameter(args[7]);
    m_n_0m = AmpParameter(args[8]);
    m_phi0_0m = AmpParameter(args[20]);
    m_theta_0m = AmpParameter(args[21]);    

  registerParameter(m_1p);
  registerParameter(m_w_1p);
  registerParameter(m_n_1p);
  registerParameter(m_phi0_1p);
  registerParameter(m_theta_1p);
  registerParameter(m_phip_1p);
  registerParameter(m_phim_1p);
  registerParameter(m_psi_1p);

  registerParameter(m_1m);
  registerParameter(m_w_1m);
  registerParameter(m_n_1m);
  registerParameter(m_phi0_1m);
  registerParameter(m_theta_1m);
  registerParameter(m_phip_1m);
  registerParameter(m_phim_1p);
  registerParameter(m_psi_1m);

  registerParameter(m_0m);
  registerParameter(m_w_0m);
  registerParameter(m_n_0m);
  registerParameter(m_phi0_0m);
  registerParameter(m_theta_0m);

  registerParameter(m_ds_ratio);

    pararray[0] = m_1p;
    pararray[1] = m_w_1p;
    pararray[2] = m_n_1p;
    pararray[3] = m_1m;
    pararray[4] = m_w_1m;
    pararray[5] = m_n_1m;
    pararray[6] = m_0m;
    pararray[7] = m_w_0m;
    pararray[8] = m_n_0m;
    pararray[9] = m_ds_ratio;
    pararray[10] = m_phi0_1p;
    pararray[11] = m_theta_1p;
    pararray[12] = m_phip_1p;
    pararray[13] = m_phim_1p;
    pararray[14] = m_psi_1p;
    pararray[15] = m_phi0_1m;
    pararray[16] = m_theta_1m;
    pararray[17] = m_phip_1m;
    pararray[18] = m_phim_1m;
    pararray[19] = m_psi_1m;
    pararray[20] = m_phi0_0m;
    pararray[21] = m_theta_0m;
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

////////////////////////////////////////////////// User Vars //////////////////////////////////
void
omegapiAngAmp::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  

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
  
       double moment[25] = {0.0};
       double calpha_array[25] = {0.0};

    for (int alpha = 0; alpha < 25; alpha++)//quantum states
	  {
	    moment[alpha] = hmoment(alpha, angvector);
	  
            calpha_array[alpha] = calpha(alpha);
        }//alpha loop
	  
  userVars[uv_Phi] = locthetaphi[2];

  userVars[uv_Pgamma] = Pgamma;
  userVars[uv_mx] = mx;

  userVars[uv_moment0] = moment[0];
  userVars[uv_moment1] = moment[1];
  userVars[uv_moment2] = moment[2];
  userVars[uv_moment3] = moment[3];
  userVars[uv_moment4] = moment[4];
  userVars[uv_moment5] = moment[5];
  userVars[uv_moment6] = moment[6];
  userVars[uv_moment7] = moment[7];
  userVars[uv_moment8] = moment[8];
  userVars[uv_moment9] = moment[9];
  userVars[uv_moment10] = moment[10];
  userVars[uv_moment11] = moment[11];
  userVars[uv_moment12] = moment[12];
  userVars[uv_moment13] = moment[13];
  userVars[uv_moment14] = moment[14];
  userVars[uv_moment15] = moment[15];
  userVars[uv_moment16] = moment[16];
  userVars[uv_moment17] = moment[17];
  userVars[uv_moment18] = moment[18];
  userVars[uv_moment19] = moment[19];
  userVars[uv_moment20] = moment[20];
  userVars[uv_moment21] = moment[21];
  userVars[uv_moment22] = moment[22];
  userVars[uv_moment23] = moment[23];
  userVars[uv_moment24] = moment[24];
  
  userVars[uv_calpha0] = calpha_array[0];
  userVars[uv_calpha1] = calpha_array[1];
  userVars[uv_calpha2] = calpha_array[2];
  userVars[uv_calpha3] = calpha_array[3];
  userVars[uv_calpha4] = calpha_array[4];
  userVars[uv_calpha5] = calpha_array[5];
  userVars[uv_calpha6] = calpha_array[6];
  userVars[uv_calpha7] = calpha_array[7];
  userVars[uv_calpha8] = calpha_array[8];
  userVars[uv_calpha9] = calpha_array[9];
  userVars[uv_calpha10] = calpha_array[10];
  userVars[uv_calpha11] = calpha_array[11];
  userVars[uv_calpha12] = calpha_array[12];
  userVars[uv_calpha13] = calpha_array[13];
  userVars[uv_calpha14] = calpha_array[14];
  userVars[uv_calpha15] = calpha_array[15];
  userVars[uv_calpha16] = calpha_array[16];
  userVars[uv_calpha17] = calpha_array[17];
  userVars[uv_calpha18] = calpha_array[18];
  userVars[uv_calpha19] = calpha_array[19];
  userVars[uv_calpha20] = calpha_array[20];
  userVars[uv_calpha21] = calpha_array[21];
  userVars[uv_calpha22] = calpha_array[22];
  userVars[uv_calpha23] = calpha_array[23];
  userVars[uv_calpha24] = calpha_array[24];
  
}

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
omegapiAngAmp::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  complex <GDouble> COne(1,0);  

   GDouble Phi = userVars[uv_Phi];
   GDouble polfrac = userVars[uv_Pgamma];
   GDouble mb1 = userVars[uv_mx];

  GDouble moment[25] = {userVars[uv_moment0], userVars[uv_moment1], userVars[uv_moment2], userVars[uv_moment3], userVars[uv_moment4], userVars[uv_moment5], userVars[uv_moment6], userVars[uv_moment7], userVars[uv_moment8], userVars[uv_moment9], userVars[uv_moment10], userVars[uv_moment11], userVars[uv_moment12], userVars[uv_moment13], userVars[uv_moment14], userVars[uv_moment15], userVars[uv_moment16], userVars[uv_moment17], userVars[uv_moment18], userVars[uv_moment19], userVars[uv_moment20], userVars[uv_moment21], userVars[uv_moment22], userVars[uv_moment23], userVars[uv_moment24]};

  GDouble calpha_array[25] = {userVars[uv_calpha0], userVars[uv_calpha1], userVars[uv_calpha2], userVars[uv_calpha3], userVars[uv_calpha4], userVars[uv_calpha5], userVars[uv_calpha6], userVars[uv_calpha7], userVars[uv_calpha8], userVars[uv_calpha9], userVars[uv_calpha10], userVars[uv_calpha11], userVars[uv_calpha12], userVars[uv_calpha13], userVars[uv_calpha14], userVars[uv_calpha15], userVars[uv_calpha16], userVars[uv_calpha17], userVars[uv_calpha18], userVars[uv_calpha19], userVars[uv_calpha20], userVars[uv_calpha21], userVars[uv_calpha22], userVars[uv_calpha23], userVars[uv_calpha24]};

 
    double wdist0 = 0.0;
    double wdist1 = 0.0;
    double wdist2 = 0.0;

          for (int alpha = 0; alpha < 25; alpha++)//quantum states
	  {
	    wdist0 += Intensity(mb1, pararray, alpha) * moment[alpha] * calpha_array[alpha];

	    wdist1 += Intensity(mb1, pararray, alpha) * moment[alpha] * calpha_array[alpha] * TMath::Sqrt(2.0) * TMath::Cos(2.0 * Phi);

	    wdist2 += Intensity(mb1, pararray, alpha) * moment[alpha] * calpha_array[alpha] * TMath::Sqrt(2.0) * TMath::Sin(2.0 * Phi);

	  }//alpha loop

	double intensity = 0.0;
	double real_sqrt_intensity = 0.0;
	double ima_sqrt_intensity = 0.0;
	
	intensity = wdist0 - polfrac * wdist1 * TMath::Cos(2.0 * Phi) - polfrac * wdist2 * TMath::Sin(2.0 * Phi);
	
       if (intensity >= 0.0) real_sqrt_intensity = sqrt(intensity);
       if (intensity < 0.0) ima_sqrt_intensity = sqrt(abs(intensity));
	
    complex <GDouble> sqrt_intensity(real_sqrt_intensity,ima_sqrt_intensity);

    complex <GDouble> amplitude = sqrt_intensity;

    //cout << real(amplitude) << "+i" << imag(amplitude) << endl;


return amplitude;
}

void omegapiAngAmp::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates  
}

#ifdef GPU_ACCELERATION
void
omegapiAngAmp::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
    
  GPUomegapiAngAmp_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
                            m_1p, m_w_1p, m_n_1p, m_1m, m_w_1m, m_n_1m,
                            m_0m, m_w_0m, m_n_0m, m_ds_ratio, m_phi0_1p,
                            m_theta_1p, m_phip_1p, m_phim_1p, m_psi_1p,
			    m_phi0_1m, m_theta_1m, m_phip_1m, m_phim_1m,
                            m_psi_1m, m_phi0_0m, m_theta_0m, useCutoff, polAngle, polFraction );

}
#endif //GPU_ACCELERATION
