//may 28th 2018, Atkinson 84, wigner D, config input
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include "UTILITIES/CobremsGeneration.hh"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "omegapiAngAmp.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "AMPTOOLS_AMPS/HelicityFrame.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

TLorentzVector targetopi(0,0,0,0.938);

double par[22];

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
  int L = lmLM[alpha][1];
  int m = lmLM[alpha][2];
  int M = lmLM[alpha][3];
  
  double normalization_num = 4.0 * TMath::Pi() * TMath::Pi();
  double normalization_denum = ((2.0 * l)+1.0) * ((2.0 * L)+1.0) * (2.0-delta(m,0)) * (2.0-delta(M,0));
  normalization = normalization_denum/normalization_num;
  
  return normalization;
}

double hmoment(int alpha, vector<double> vector)
{
  double loccostheta1 = TMath::Cos(vector[0]);
  double locsintheta1 = TMath::Sin(vector[0]);
  double locphi1 = vector[1];
  double loccosthetaH1 = TMath::Cos(vector[2]);
  double locsinthetaH1 = TMath::Sin(vector[2]);
  double locphiH1 = vector[3];

  int l = lmLM[alpha][0];
  int m = lmLM[alpha][1];
  int L = lmLM[alpha][2];
  int M = lmLM[alpha][3];

   double moment = 0.0;
  if (alpha < 15)
  {
    moment = 0.5 * std::real(wignerD(L, M, m, loccostheta1, locphi1) * wignerD(l, m, 0, loccosthetaH1, locphiH1) + pow(-1.0,L+M) * wignerD(L, M, m, loccostheta1, locphi1) * wignerD(l, m, 0, loccosthetaH1, locphiH1));
  }
  else {moment = 0.5 * std::real(wignerD(L, M, m, loccostheta1, locphi1) * wignerD(l, m, 0, loccosthetaH1, locphiH1) - pow(-1.0,L+M) * wignerD(L, M, m, loccostheta1, locphi1) * wignerD(l, m, 0, loccosthetaH1, locphiH1));}

  return moment;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
omegapiAngAmp::omegapiAngAmp( const vector< string >& args ):
  UserAmplitude< omegapiAngAmp >( args ),
  m_fastCalc( false )
{

    assert( args.size() == 23);

    polAngle = AmpParameter( args[0] );

    registerParameter( polAngle );


    ds_ratio = atof( args[10].c_str() );
    
    //1+ state
    m_1p = atof( args[1].c_str() );
    w_1p = atof( args[2].c_str() );
    n_1p = atof( args[3].c_str() );
    phi0_1p = atof( args[11].c_str() );
    theta_1p = atof( args[12].c_str() );
    phip_1p = atof( args[13].c_str() );
    phim_1p = atof( args[14].c_str() );
    psi_1p = atof( args[15].c_str() );
    //1- state    
    m_1m = atof( args[4].c_str() );
    w_1m = atof( args[5].c_str() );
    n_1m = atof( args[6].c_str() );
    phi0_1m = atof( args[16].c_str() );
    theta_1m = atof( args[17].c_str() );
    phip_1m = atof( args[18].c_str() );
    phim_1m = atof( args[19].c_str() );
    psi_1m = atof( args[20].c_str() );
    //0- state    
    m_0m = atof( args[7].c_str() );
    w_0m = atof( args[8].c_str() );
    n_0m = atof( args[9].c_str() );
    phi0_0m = atof( args[21].c_str() );
    theta_0m = atof( args[22].c_str() );
    

    par[0] = m_1p;
    par[1] = w_1p;
    par[2] = n_1p;
    par[3] = m_1m;
    par[4] = w_1m;
    par[5] = n_1m;
    par[6] = m_0m;
    par[7] = w_0m;
    par[8] = n_0m;
    par[9] = ds_ratio;
    par[10] = phi0_1p;
    par[11] = theta_1p;
    par[12] = phip_1p;
    par[13] = phim_1p;
    par[14] = psi_1p;
    par[15] = phi0_1m;
    par[16] = theta_1m;
    par[17] = phip_1m;
    par[18] = phim_1m;
    par[19] = psi_1m;
    par[20] = phi0_0m;
    par[21] = theta_0m;
    
  mG0_omega = 0.0085;
  mG0_b1 = 0.143;
    // Initialize coherent brem table
    // Do this over the full range since we will be using this as a lookup
    float Emax  = 12.0;
    float Epeak = 9.0;
    float Elow  = 0.139*2;
    float Ehigh = 12.0;

    int doPolFlux=0;  // want total flux (1 for polarized flux)
    float emitmr=10.e-9; // electron beam emittance
    float radt=50.e-6; // radiator thickness in m
    float collDiam=0.005; // meters
    float Dist = 76.0; // meters
    CobremsGeneration cobrems(Emax, Epeak);
    cobrems.setBeamEmittance(emitmr);
    cobrems.setTargetThickness(radt);
    cobrems.setCollimatorDistance(Dist);
    cobrems.setCollimatorDiameter(collDiam);
    cobrems.setCollimatedFlag(true);
    cobrems.setPolarizedFlag(doPolFlux);

    // Create histogram
    totalFlux_vs_E = new TH1D("totalFlux_vs_E", "Total Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
    polFlux_vs_E   = new TH1D("polFlux_vs_E", "Polarized Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
    polFrac_vs_E   = new TH1D("polFrac_vs_E", "Polarization Fraction vs. E_{#gamma}", 1000, Elow, Ehigh);

    // Fill totalFlux
    for(int i=1;i<=totalFlux_vs_E->GetNbinsX(); i++){
        double x = totalFlux_vs_E->GetBinCenter(i)/Emax;
        double y = 0;
        //if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
        y = cobrems.Rate_dNtdx(x);
        totalFlux_vs_E->SetBinContent(i, y);
    }

    doPolFlux=1;
    cobrems.setPolarizedFlag(doPolFlux);
    // Fill totalFlux
    for(int i=1;i<=polFlux_vs_E->GetNbinsX(); i++){
        double x = polFlux_vs_E->GetBinCenter(i)/Emax;
        double y = 0;
        //if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
        y = cobrems.Rate_dNcdx(x);
        polFlux_vs_E->SetBinContent(i, y);
    }

    polFrac_vs_E->Divide(polFlux_vs_E, totalFlux_vs_E);
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
  return sumH;
}

//Define t_star (sum of production density matrix)
std::complex<double> t_star_LM(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int ialpha, int ibeta, int L, int M) {
  return TMath::Sqrt((2.0*J_spin(ibeta) + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * HelicitySum(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, ialpha, ibeta, L, M);
}

double SingleIntensity(double x, double *par, int l, int m, int L, int M, int ialpha, int ibeta) {
    double phiplus_alpha = 0;
    double phiminus_alpha = 0;
    double psi_alpha = 0;
    double phiplus_beta = 0;
    double phiminus_beta = 0;
    double psi_beta = 0;
    if (ialpha < 2) {
      phiplus_alpha = par[5*ialpha + 12];
      phiminus_alpha = par[5*ialpha + 13];
      psi_alpha = par[5*ialpha + 14];
    }
    if (ibeta < 2) {
      phiplus_beta = par[5*ibeta + 12];
      phiminus_beta = par[5*ibeta + 13];
      psi_beta = par[5*ibeta + 14];
    }
  return std::real(t_star_LM(par[5*ialpha + 10], par[5*ialpha + 11], phiplus_alpha, phiminus_alpha, psi_alpha, par[5*ibeta + 10], par[5*ibeta + 11], phiplus_beta, phiminus_beta, psi_beta, ialpha, ibeta, L, M) * f_Llm(x, par[3*ialpha + 0], par[3*ialpha + 1], par[3*ialpha + 2], ialpha, par[3*ibeta + 0], par[3*ibeta + 1], par[3*ibeta + 2], ibeta, par[9], l, m, L) * clebschGordan(1, l, 0, 0, 1, 0));
}

double Intensity(double x, double *par, int i){
  
    double IntensitySum = 0;
  for (int ialpha = 0; ialpha < 3; ialpha++) { //double sum over combination of JP states
    for (int ibeta = 0; ibeta < 3; ibeta++) {
      if ( (ialpha != ibeta) & (eta_parity(ialpha) == eta_parity(ibeta)) ) //Only want interference moments with opposite parity
	continue;
      IntensitySum += SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
    }
  }
  return IntensitySum;
}

double wdist(int k, double x, double Phi, vector<double> vector)
  {
    double dist = 0.0;

          for (int alpha = 0; alpha < 25; alpha++)//quantum states
	  {
	    if (k == 0){dist += Intensity(x, par, alpha) * hmoment(alpha, vector) * calpha(alpha);}

	    if (k == 1){dist += Intensity(x, par, alpha) * TMath::Sqrt(2.0) * TMath::Cos(2.0 * Phi) * hmoment(alpha, vector) * calpha(alpha);}

	    if (k == 2){dist += Intensity(x, par, alpha) * TMath::Sqrt(2.0) * TMath::Sin(2.0 * Phi) * hmoment(alpha, vector) * calpha(alpha);}

	  }//alpha loop
    //cout << "dist = " << dist << endl;

    return dist;
  }

double sqrtIntensity(double polfrac, double x, double Phi, vector<double> vector)
  {
    double intensity = 0.0;

	intensity = wdist(0, x, Phi,vector) + polfrac * wdist(1, x, Phi,vector) * TMath::Cos(2.0 * Phi) +  polfrac * wdist(2, x, Phi,vector) * TMath::Sin(2.0 * Phi);
	
    //cout << "sqrtintensity = " << TMath::Sqrt(intensity) << endl;
    return TMath::Sqrt(intensity);
  }

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
omegapiAngAmp::calcAmplitude( GDouble** pKin ) const
{
  bool useCutoff=true;
  complex <GDouble> i(0, 1), COne(1, 0),CZero(0,0);  

  ////////////////////////////////////////////////// Get Lorentz Vectors of Particles //////////////////////////////////

  TLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  TLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

  GDouble m0_omega=0.783;


  //Exprected particle list: 
  // b1(pi0 omega(pi0 "rho"(pi- pi+)))
  //    2         3         4   5

  TLorentzVector rhos_pim(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector rhos_pip(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
  TLorentzVector rho = rhos_pip + rhos_pim;

  if( useCutoff && rho.M()+0.135 > m0_omega+3.0*mG0_omega){
    ////cout << "s";
    return CZero;
  }

  TLorentzVector omegas_pi(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector omega = rho + omegas_pi;

  if(useCutoff && fabs(omega.M()-m0_omega) > 3.0*mG0_omega){
    ////cout << "s";
    return CZero; 
  }

  TLorentzVector Xs_pi(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);
  TLorentzVector X = omega + Xs_pi;

    ////////////////////////////////////////////////// Boost Particles and Get Angles//////////////////////////////////

  // orientation of production plane in lab
  //GDouble alpha = recoil.Vect().Phi();
  
  //Helicity coordinate system
  TLorentzVector Gammap = beam + targetopi;
 
  vector <double> locthetaphi = getthetaphi(omega, X, beam, Gammap);
  
  vector <double> locthetaphih = getthetaphi(rhos_pip, omega, X, X, rhos_pim, Gammap); 

  vector <double> angvector{locthetaphi[0], locthetaphi[1], locthetaphih[0], locthetaphih[1]};
  
  TLorentzVector z_rf = X;
  TLorentzVector InverseOfX_rf = beam;
  TVector3 zrfboost = Gammap.BoostVector();
// boost x to zrf
  InverseOfX_rf.Boost(-1.0*zrfboost);
  //boost z to rf
  z_rf.Boost(-1.0*zrfboost);

  //particle 1 (positively charged)
  TLorentzVector particle1_rf = omega;
  //boost particle 1 to zrf
  particle1_rf.Boost(-1.0*zrfboost);
  //boost particle 1 to z
  TVector3 zboost = z_rf.BoostVector();
  TLorentzVector particle1_z = particle1_rf;
  particle1_z.Boost(-1.0*zboost);
  TVector3 particle_z = particle1_z.Vect();

  //get the unit vectors in space
  TVector3 particle_zunit = (particle_z).Unit();
  TVector3 z_rfunit = (z_rf.Vect()).Unit();
  TVector3 InverseOfX_rfunit = (InverseOfX_rf.Vect()).Unit();

  //calculate phi
  TVector3 y = ((InverseOfX_rfunit).Cross(z_rfunit)).Unit();
  TVector3 x = (y.Cross(z_rfunit)).Unit();
  TVector3 z = z_rfunit;

    TVector3 eps(cos(polAngle), sin(polAngle), 0.0); // beam polarization vector
    GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

    // vector meson production from K. Schilling et. al.
    GDouble Pgamma;
    if(polFraction >= 0.) Pgamma = polFraction;
    else{
       int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
       if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
          Pgamma = 0.;
       }
       else Pgamma = polFrac_vs_E->GetBinContent(bin);
    }


   double mx = X.M();
    complex <GDouble> amplitude = sqrtIntensity(Pgamma, mx, Phi, angvector);

   
  return amplitude;

}
