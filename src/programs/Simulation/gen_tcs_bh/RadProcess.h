#ifndef RadProcess_h
#define RadProcess_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TF1.h"
//#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Options_tcs.h"
#include "TreatOptions.h"
#include "Utils.h"
#include <unistd.h>
#include "Constants.h"
#include "TableProvider.h"
#include "PDGinfo.h"

class RadProcess {

   public: 

	float EqRad_lenght(float);
	float SampleNph(float Nph_steps[],float E, float AA, float ZZ, float d);
	TLorentzVector ElectronRC(TLorentzVector el_in, float Ecut, float Z, float d, float tot, int &rad) ;
	float IntegralBr(float Eel, float Emin, float Emax,float AA, float ZZ, float LL); 
	float IntegralNph(float E, float kmin, float kmax,  float d);
	double Thsample(double Ei);
	float Eg_rdint(float Eel,float QQ);
	float Eg_brem(float Eel, float ecut,float LL);
	double VFluxFactor(double Q2_max, double E, double nu);
	double BremstrahlungSpectraNgammadiff(double E,double nu,double ZZ,double AA, double dtarget);
	double thetag_brem(double Ei, double Ef, double ZZ);
	float InitBrmProfile(float Eel, float Emin, float Emax); 
	float BeamProfileRescale_bmr(float Egam, float Eel, float ib, float ,float); 
	TLorentzVector TLBrDeviate(TLorentzVector el_in, float Z, float Ecut, float LL);
	float UNX (double,double);

   private:

};

#endif



