#ifndef Kinematics_h
#define Kinematics_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <unistd.h>
#include "Options_tcs.h"
class Kinematics {

	public:
		float fCosTh(float, float,float, float, float,float, float);
		float ftmin(float, float, float, float);
		float fmass(float,float);
		float fen(float,float);
		float fmom(float,float);
		float fEg_out_lab(float, float, float);
		float fQ2 (float, float,float);
		float fXbj(float, float);
		float fXbj(float, float,float);
		float Wsqr(float, float,float);
	private:

};

#endif



