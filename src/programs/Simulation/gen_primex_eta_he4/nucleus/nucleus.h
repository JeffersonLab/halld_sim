#include "AV18_deut.h"

#include "He_p.h"
#include "He_n.h"

#include "C12_MF_New.hh"

//#include "cross_sections.h"
#include "masses.h"

const double GeVfm = 0.1973;

class nucleus {
  
 private:
  TRandom3 * myRand;
  
  double kF = 0.250;

  double sq(double x)
  { return x*x; }

  double get_phiSq(double *phiPtr, double k_rel)
  {

    double bin = k_rel / GeVfm / 0.1;

    if (bin < 0.)
      return 0.;
    if (bin < 1.)
      return bin * phiPtr[0];
    if (bin > 100.)
      return 0.;

    int b = bin;
    double x = bin - b;
    return (x*phiPtr[b] + (1.-x)*phiPtr[b-1]) / pow(GeVfm,3);

  }

  double get_SF_He_p(double k, double E)
  {
    
    double binx = k / GeVfm / 0.05 + 0.5;
    double biny = E * 1000 + 0.5;
    
    if (k < 0.)
      return 0.;
    if (E < 0.)
      return 0.;
    if (binx < 1.)
      binx = 1.;
    if (biny < 1.)
      biny = 1.;
    if (binx > 200.)
      return 0.;
    if (binx > 1000.)
      return 0.;
    
    int bx = binx;
    int by = biny;
    
    double x = binx - bx;
    double y = biny - by;
    
    double f1 = x*SF_He_p_array[bx][by-1] + (1.-x)*SF_He_p_array[bx-1][by-1];
    double f2 = x*SF_He_p_array[bx][by] + (1.-x)*SF_He_p_array[bx-1][by];
    
    return (y*f2 + (1.-y)*f1)*1000/pow(2*M_PI*GeVfm,3);
    
  }
  
  double get_SF_He_n(double k, double E)
  {
    
    double binx = k / GeVfm / 0.05 + 0.5;
    double biny = E * 1000 + 0.5;
    
    if (k < 0.)
      return 0.;
    if (E < 0.)
      return 0.;
    if (binx < 1.)
      binx = 1.;
    if (biny < 1.)
      biny = 1.;
    if (binx > 200.)
      return 0.;
    if (binx > 1000.)
      return 0.;
    
    int bx = binx;
    int by = biny;
    
    double x = binx - bx;
    double y = biny - by;
    
    double f1 = x*SF_He_n_array[bx][by-1] + (1.-x)*SF_He_n_array[bx-1][by-1];
    double f2 = x*SF_He_n_array[bx][by] + (1.-x)*SF_He_n_array[bx-1][by];
    
    return (y*f2 + (1.-y)*f1)*1000/pow(2*M_PI*GeVfm,3);
    
  }

  double get_SF_C12(double k, double E)
  {
  
    double binx = k / GeVfm / 0.03 + 0.5;
    double biny = E * 1000 / 0.5 + 0.5;

    if (k < 0.)
      return 0.;
    if (E < 0.)
      return 0.;
    if (binx < 1.)
      binx = 1.;
    if (biny < 1.)
      biny = 1.;
    if (binx > 200.)
      return 0.;
    if (binx > 1000.)
      return 0.;
  
    int bx = binx;
    int by = biny;

    double x = binx - bx;
    double y = biny - by;

    double f1 = x*C12_MF_New[bx][by-1] + (1.-x)*C12_MF_New[bx-1][by-1];
    double f2 = x*C12_MF_New[bx][by] + (1.-x)*C12_MF_New[bx-1][by];

    return (y*f2 + (1.-y)*f1)*1000/pow(2*M_PI*GeVfm,3);
  
  }

 public:
  nucleus(TRandom3 * thisRand)
    { myRand = thisRand; }

  ~nucleus()
    { }

  void sample_SF_deut_p(double &weight, double &k, double &E)
  {

    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);
    
    // E is fixed for deuterium
    E = mP - (m_2H - sqrt(mN*mN + k*k));

    // Calling spectral function
    weight *= get_phiSq(AV18_deut,k)/pow(2.*M_PI,3);

  }

  void sample_SF_deut_n(double &weight, double &k, double &E)
  {

    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);
    
    // E is fixed for deuterium
    E = mN - (m_2H - sqrt(mP*mP + k*k));

    // Calling spectral function
    weight *= get_phiSq(AV18_deut,k)/pow(2.*M_PI,3);

  }
  
  void sample_SF_He_p(double &weight, double &k, double &E)
  {
    
    //cout << "weight0 " << weight << endl;
    // Random sampling of k
    k = myRand->Exp(kF);
    //cout << "kF " << kF << " k " << k << endl;
    weight *= kF*exp(k/kF);
    //cout << "weight1 " << weight << endl;
    weight *= 4*M_PI*sq(k);
    //cout << "weight2 " << weight << endl;
    // Random sampling of E
    double Emin = sq(k)/(2.*m_3H);
    double Emax = Emin + 0.04;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);
    //cout << "weight3 " << weight << endl;
    // Calling spectral function
    weight *= get_SF_He_p(k,E);
    //cout << "weight4 " << weight << endl;
  }

  void sample_SF_He_d(double &weight, double &k, double &E)
  {
    
    //cout << "weight0 " << weight << endl;
    // Random sampling of k
    k = myRand->Exp(kF);
    //cout << "kF " << kF << " k " << k << endl;
    weight *= kF*exp(k/kF);
    //cout << "weight1 " << weight << endl;
    weight *= 4*M_PI*sq(k);
    //cout << "weight2 " << weight << endl;
    // Random sampling of E
    double Emin = sq(k)/(2.*m_2H);
    double Emax = Emin + 0.04;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);
    //cout << "weight3 " << weight << endl;
    // Calling spectral function
    weight *= get_SF_He_p(k,E);
    //cout << "weight4 " << weight << endl;
  }
  

  void sample_SF_He_n(double &weight, double &k, double &E)
  {
    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);

    // Random sampling of E
    double Emin = sq(k)/(2.*m_3He);
    double Emax = Emin + 0.04;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);

    // Calling spectral function
    weight *= get_SF_He_n(k,E);
    //cout <<"in nucleus n weight " << weight << " k " << k << " E " << E << endl;
  }

  void sample_SF_C12_p(double &weight, double &k, double &E)
  {
    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);

    // Random sampling of E
    double Emin = sq(k)/(2.*m_11B);
    double Emax = Emin + 0.3;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);

    // Calling spectral function
    weight *= get_SF_C12(k,E);

  }


  void sample_SF_C12_n(double &weight, double &k, double &E)
  {
    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);

    // Random sampling of E
    double Emin = sq(k)/(2.*m_11C);
    double Emax = Emin + 0.3;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);

    // Calling spectral function
    weight *= get_SF_C12(k,E);

  }
  
  void sample_SF_C12_d(double &weight, double &k, double &E)
  {
    // Random sampling of k
    k = myRand->Exp(kF);
    weight *= kF*exp(k/kF);
    weight *= 4*M_PI*sq(k);
    
    // Random sampling of E
    double Emin = sq(k)/(2.*m_2H);
    double Emax = Emin + 0.3;
    E = Emin + (Emax - Emin)*myRand->Uniform();
    weight *= (Emax - Emin);

    // Calling spectral function
    weight *= get_SF_C12(k,E);

  }

  
};
