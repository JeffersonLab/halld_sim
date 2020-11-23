#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>
#include <complex>
#include <string>
#include <time.h>

#include "IUAmpTools/FitResults.h"
#include "TFile.h"

//#include "wave.h"
//#include "3j.h"

#include "TFile.h"

using namespace std;




// UNITS = GEV
#define MP      0.93827203
#define MPI     0.13957061
#define META    0.547682
#define METAP   0.95778
#if !defined( M_PI )
#define M_PI 3.141592654    
#endif
#define LMAX    2             // highest wave





std::complex<double> pw_refl(string fitresFile, int L, int hel[3]);
std::complex<double> sdme_refl(int alp, int L1, int L2, int M1, int M2, string fitFile);
std::complex<double> Moments_refl(int alp, int L, int M, string file_fit);
double clebsch(double j1, double j2, double j3, double m1, double m2);








int main( int argc, char* argv[] ){
    
    // these params should probably come in on the command line
    double lowMass = 0.7;
    double highMass = 3.0;
    enum{ kNumBins = 45 };
    double lowt = 0;
    double hight = 1.2;
    enum{ kNumBinst = 4 };
    string fitDir("etaprimepi0");
                     
    // set default parameters
    
    string outfileName("");
    
    // parse command line
    
    for (int i = 1; i < argc; i++){
        
        string arg(argv[i]);
        
        if (arg == "-o"){
            if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
            else  outfileName = argv[++i]; }
        if (arg == "-h"){
            cout << endl << " Usage for: " << argv[0] << endl << endl;
            cout << "\t -o <file>\t Ouput text file" << endl;
            exit(1);}
    }
    
    if (outfileName.size() == 0){
        cout << "No output file specified" << endl;
        exit(1);
    }
    


 


  double step = ( highMass - lowMass ) / kNumBins;
    double stept = ( hight - lowt ) / kNumBinst;
    
    ofstream outfile;
    outfile.open( outfileName.c_str() );

    outfile <<"M"<<"\t"<<"t"; //First line contains names of variables, first two colomns correspond to M(invariant mass) and t

    for (int L = 0; L<= pow(LMAX,2); L++) {// the rest of the colomn correspond to moments
      for (int M = 0; M<= L; M++) {
    
	outfile<<"\t"<<"H0_"<<L<<M<<"\t"<<"H0_"<<L<<M<<"uncert.";
	outfile<<"\t"<<"H1_"<<L<<M<<"\t"<<"H1_"<<L<<M<<"uncert.";
      }}
    outfile<<endl;



    // descend into the directory that contains the bins
    chdir( fitDir.c_str() );



    //Looping through M_eta_pi and t bins
    for( int i = 0; i < kNumBins;i++ ){
      for( int j = 0; j < kNumBinst; j++ ){  
	cout<<"bin "<<i<<"_"<<j<<endl;
        
        ostringstream dir;
        dir << "bin_" << i<<"_"<<j;
        chdir( dir.str().c_str() );
        
        string resultsFile;
        resultsFile ="bin_" + std::to_string(i)+"_"+ std::to_string(j)+".fit";

        
        // print out the bin center
        outfile << lowMass + step * i + step / 2. << "\t";
	outfile << lowt + stept * j + stept / 2. << "\t";

  
    for (int L = 0; L<= pow(LMAX,2); L++) {// calculating moments and writing to a file
    for (int M = 0; M<= L; M++) {
     
      outfile << real(Moments_refl(0, L, M, resultsFile))<< "\t"<< 0 <<"\t";
      outfile << real(Moments_refl(1, L, M, resultsFile))<< "\t"<< 0 <<"\t";
      // if(L==0 && M==0 && i==44 && j==3)cout <<" H0_00=  "<< real(Moments_refl(0, L, M, resultsFile))<<" H1_00=  "<< real(Moments_refl(1, L, M, resultsFile))<<endl;
        }}
        outfile << endl;
        
        chdir( ".." );
    }
    }




    return 0;
}











std::complex<double> pw_refl(string fitresFile, int L, int hel[3]){
  /* returns partial waves for g p --> (eta pi)_L p
  * partial waves are in the reflectivity basis
  * fitresFile is the name of the fit output file for given m_eta_pi and t bin
  * hel = {epsilon, m, k}
  * epsilon is the relfectivity
  * k = 0,1 is the nucleon non-flip (0) or flip (1)
  */


 std::complex<double> pw (0.0,0.0), ui (0.,1.), zero(0.0, 0.0);
  
  int eps = hel[0], m = hel[1], k = hel[2];


FitResults fitres(fitresFile);

 if( !fitres.valid() )return zero;
 // eps = +/- 1 ; k = 0,1 ; |m|<= L
  if( abs(eps)!=1 || abs(2*k-1)!=1 || abs(m)>L) {
    return zero;
  }
  
  
  // only positive refelctivity
  if (eps == -1){
    return zero;
  }
  // only nucleon non-flip
  if (k == 1){
    return zero;
  }
  // only positive m
  if (m < 0){
    return zero;
  }
 
  if (eps == 1){
    switch (L) {
      case 0:
        // a0(980)
	pw=fitres.scaledProductionParameter("EtaPrimePi0::PositiveRe::S0+");
        break;
      case 1:
        // pi1(1600)
	if(m==0)pw=fitres.productionParameter("EtaPrimePi0::PositiveRe::P0+");
	else pw=fitres.productionParameter("EtaPrimePi0::PositiveRe::P1+");
	break;
      case 2:
        // a2(1320) + a2(1700)
        if(m==0)pw=fitres.productionParameter("EtaPrimePi0::PositiveRe::D0+");
	else if(m==1)pw=fitres.productionParameter("EtaPrimePi0::PositiveRe::D1+");
	else pw=fitres.productionParameter("EtaPrimePi0::PositiveRe::D2+");
	break;
        
      default:
        pw = zero;
        break;
    }
  
  
  }
  return pw;
}






std::complex<double> sdme_refl(int alp, int L1, int L2, int M1, int M2, string fitFile){ // code formula (D8) from 10.1103/PhysRevD.100.054017 //alp corresponds to H0 or H1
  // fitresFile contains fit results in given {t,m_etapi} bin
  std::complex<double> rho (0.0,0.0), ui (0.0, 1.0);
  std::complex<double> pw1 (0.0,0.0), pw2 (0.0,0.0);
  int hel1[3], hel2[3]; // hel = [eps, m , k]
  std::complex<double> fac (1.0,0.0);
  
  // sum over relfectivity and nucleon (non-)flip
  for (int e = -1; e <= 1; e +=2) {
  for (int k =  0; k <= 1;   k++) {
    hel1[0] = e; hel2[0] = e;
    hel1[2] = k; hel2[2] = k;
    
    switch (alp) {   //alp corresponds to H0 or H1
      case 0:
        fac = 1.0;
        hel1[1] = M1; hel2[1] = M2;
        pw1 = pw_refl(fitFile, L1, hel1); //  fitresFile contains fitresults for given bin in { t_{pp}, m_{eta pi}}, L , hel = {epsilon, m proj. of L, k}
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        
        fac = pow(-1.,M1-M2);
        hel1[1] = -M1; hel2[1] = -M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
       
        break;
        
      case 1:
        fac = -(double)e*pow(-1.,M1);
        
        hel1[1] = -M1; hel2[1] = M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        
        fac = -(double)e*pow(-1.,M2);
        hel1[1] = M1; hel2[1] = -M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        break;
        
      case 2:
        fac = -ui*(double)e*pow(-1.,M1);
        
        hel1[1] = -M1; hel2[1] = M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        
        fac = +ui*(double)e*pow(-1.,M2);
        hel1[1] = M1; hel2[1] = -M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        break;
        
      case 3:
        fac = 1.0;
        
        hel1[1] = M1; hel2[1] = M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        
        fac = -pow(-1.,M1-M2);
        hel1[1] = -M1; hel2[1] = -M2;
        pw1 = pw_refl(fitFile, L1, hel1);
        pw2 = pw_refl(fitFile, L2, hel2);
        rho += fac * pw1 * conj(pw2);
        break;
        
      default:
        return 0.0; break;
    }
  }}
  // no phase space added
  return rho;
}





std::complex<double> Moments_refl(int alp, int L, int M, string file_fit){ //H0 or H1 , L, M, fitresults file for given t and  invariant mass bin 
  // WARNING: the sum extends to max(l1,l2) = LMAX 
  std::complex<double> mom = 0.0;
  std::complex<double> rho = 0.0;
  double cg1, cg2;
  int fac = -1; if(alp==0){fac = +1;};

  // sum over angular momenta l1 and l2     //based on A9 parity relation
  // and m1 and m2 (m1=M+m2)
  for (int l1 = 0;l1 <= LMAX ; l1 +=1 ) {
  for (int l2 = 0; l2 <= LMAX ; l2 +=1 ) {
      for (int m2 = -l2; m2 <= l2; m2 +=1 ) {
        cg1 = clebsch(l2,L,l1,0,0);   // m1,m2 and M are =0
        cg2 = clebsch(l2,L,l1,m2,M);  // 6th argument m1=M+m2
        rho = sdme_refl(alp, l1, l2, m2+M,m2, file_fit);
        mom += fac*sqrt( (2.0*l2+1)/(2.0*l1+1) )*cg1*cg2*rho;
	
      }}}
 
  return mom;
}






int isfrac(double x){
  int val = 1;
  double intpart;
  if( modf(x, &intpart) == 0.0){ val = 0;}
  return val;
}


long int factorial(int n)
{
  if(n > 1)
    return n * factorial(n - 1);
  else
    return 1;
}




double clebsch(double j1, double j2, double j3, double m1, double m2){
  double cg = 0.0;
  double m3 = m1+m2;
  
  // check input
  if (isfrac(j1+j2+j3) || isfrac(j1+m1) || isfrac(j2+m2) || isfrac(j3+m3) ||
      isfrac(j3-j1-m2) || isfrac(j3-j2+m1)) {
    return 0.0;
  }
  
  
  // Check for conditions that give CG = 0.
  if ( j3 < fabs(j1-j2) || j3 < fabs(j1-j2)||
      fabs(m1) > j1 || fabs(m2) > j2 || fabs(m3) > j3 ) {
    return 0.0;
  }
  
  // Compute the Clebsch-Gordan coefficient
  cg = sqrt( (2*j3+1)/ factorial( round(j1+j2+j3+1) ) );
  cg = cg * sqrt( factorial(round(j1+j2-j3))*factorial(round(j2+j3-j1))*factorial(round(j3+j1-j2))   );
  cg = cg * sqrt( factorial(round(j1+m1))*factorial(round(j1-m1)) );
  cg = cg * sqrt( factorial(round(j2+m2))*factorial(round(j2-m2)) );
  cg = cg * sqrt( factorial(round(j3+m3))*factorial(round(j3-m3)) );
  
  double sum = 0.0, term;
  for (int k = 0; k < 99; k++) {
    if (j1+j2-j3-k < 0) continue;
    if (j3-j1-m2+k < 0) continue;
    if (j3-j2+m1+k < 0) continue;
    if (j1-m1-k    < 0) continue;
    if (j2+m2-k    < 0) continue;
    term = factorial(round(j1+j2-j3-k)) * factorial(round(j3-j1-m2+k)) ;
    term = term * factorial(round(j3-j2+m1+k)) * factorial(round(j1-m1-k));
    term = term * factorial(round(j2+m2-k)) * factorial(k);
    term = term * pow(-1, k);
    sum = sum + 1./term;
  }
  
  cg = cg * sum;
  
  return cg;
}




