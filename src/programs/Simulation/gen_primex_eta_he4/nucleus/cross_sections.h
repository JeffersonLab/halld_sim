#include <iostream>
#include <cmath>
#include "masses.h"

using namespace std;

class cross_sections {
  
 private:
  
  double sq(double x)
  {
    return x*x;
  }
  
  double dipole_F(double t, double tmin, double tmax)
  {
    //double Lambda = 0.71;
    double Lambda = 1.346185634837375;
    
    return pow(sq(Lambda) - t,-4.) * 3./(pow(sq(Lambda) - tmin,-3.) - pow(sq(Lambda) - tmax,-3.));
  }
  
 public:
  cross_sections() {
  }

  ~cross_sections() {
  }

  double sigma_meson_simple(double s, double cosThetaCM)
  {
    return pow(1.02-cosThetaCM,-5)*pow(1.02+cosThetaCM,-4)*pow(s,-7);
  }
  
  double sigma_pi0_p(double s, double t)
  {
    double a = 14.79;
    double b = 7.37;
    double c = 7.49;
    double d = 5.36;
    
    double m_m = mpi0;
    double m_B = mN;
    
    double A = s;
    double B = sq(m_B) - sq(m_m) - s;
    double C = sq(m_m);
    double D = sq(B) - 4*A*C;
    
    double k = (s - sq(mN))/(2.*mN);
    double tmax = sq(m_m) - k*mN/s*(- B - sqrt(D));
    double tmin = sq(m_m) - k*mN/s*(- B + sqrt(D));
    double x = (tmax - t)/(tmax - tmin);
    
    return pow(a/s,b)*pow(x,c)*exp(d*sq(log(x)));
  }
  
  double sigma_pip_n(double s, double cosThetaCM)
  {
    
    double A = 9.490;
    double b = 5.329;
    double c = 4.638;
    
    return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
  }
  
  double sigma_pim_p(double s, double cosThetaCM)
  {
    
    double A = 10.240;
    double b = 5.329;
    double c = 4.638;
    
    return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
  }
  
  
  double sigma_rho0_p_old(double s, double cosThetaCM)
  {
    const double b=-3.7;
    const double c=-2.2;
    const double a=5.82005e7;
    
    return 0.75*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
  }
  
  double sigma_rho0_p(double s, double t, double cosThetaCM)
  {
    
    const double A = 9.00572806e+04;
    const double B = 5.94050683e+00;
    const double C = 1.09653763e+08;
    const double D = 3.83373851e+00;
    const double E = 1.79523684e+00;
    
    return A*exp(B*t)+C*pow(s,-7.)*pow(1.2-cosThetaCM,-D)*pow(1.05+cosThetaCM,-E);
  }
  
  double sigma_rhom_p(double s, double t, double cosThetaCM)
  {
    
    //return (1.1/15.9)*sigma_rho0_p(s,t,cosThetaCM);                                                                                                                                                                                                                             
    //return 1;                                                                                                                                                                                                                                                                   
    return pow(1-cosThetaCM,-3)*pow(s,-7);
  }
  
  double sigma_omega_p_old(double s, double cosThetaCM)
  {
    const double b=-3.7;
    const double c=-2.2;
    const double a=5.82005e7;
    
    return 0.25*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
  }
  
  double sigma_omega_p(double s, double t, double cosThetaCM)
  {
    
    return sigma_rho0_p(s,t,cosThetaCM)/3.;
  }
  
  double sigma_phi_p(double s, double t, double cosThetaCM)
  {
    
    const double A = 977.2505868903087;
    const double B = 3.069753904605845;
    const double C = 1.561280750635443;
    const double D = 0.16493020208550144;
    const double E = 1320143.8506911423;
    
    return A*exp(B*t)+E*pow(s,-7.)*exp(C*sq(cosThetaCM-D));
  }
  
  double sigma_phi_n(double s, double t, double cosThetaCM)
  {
    
    return sigma_phi_p(s, t, cosThetaCM);
  }
  
  double sigma_Jpsi_p(double s, double t)
  {
    //double sig0 = 11.3; //nb
    //double beta = 1.3;

    double sig0 = 5.850224703362592;
    double beta = 1.1864734541248247;
    
    double x = (sq(mJpsi) + 2*mN*mJpsi)/(s - sq(mN));
    if ((x > 1.) or (x < 0.))
      return 0.;
    
    double sig = sig0 * pow(1. - x,beta);
    
    double pi = (sq(mN) - s)/(2.*sqrt(s));
    double pf = 0.5*sqrt(sq(sq(mN) - sq(mJpsi))/s - 2.*(sq(mN) + sq(mJpsi)) + s);
    double tmin = sq(sq(mJpsi))/(4.*s) - sq(pi + pf);
    double tmax = sq(sq(mJpsi))/(4.*s) - sq(pi - pf);
    
    return sig * dipole_F(t, tmin, tmax);
  }
  
  double sigma_Jpsi_p(double s, double t, double QSq)
  {
    double n = 2.44;
    
    return pow(sq(mJpsi)/(QSq + sq(mJpsi)),n) * sigma_Jpsi_p(s, t);
  }
  
  double R_Jpsi_p(double QSq)
  {
    double a = 2.164;
    double n = 2.131;
    return pow((a*sq(mJpsi) + QSq)/(a*sq(mJpsi)),n) - 1.;
  }
  
  double sigma_deltapp_pim(double s, double cosThetaCM)
  {
    double a =50955200;
    double b =2.93657;
    double c =1.74578;
    
    return pow(s,-7)*a*pow(1-cosThetaCM,-b)*pow(1+cosThetaCM,-c);
  }
  
  double sigma_deltap_pim(double s, double cosThetaCM)
  {
    
    return sigma_deltapp_pim(s,cosThetaCM);
  }

};
