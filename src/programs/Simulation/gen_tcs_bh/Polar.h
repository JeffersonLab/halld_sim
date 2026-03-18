#ifndef Polar_h
#define Polar_h
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
class Polar {

	public:
		float cross_poldilut(float p, float m, float spindir, float dil);
		float cross_doublepoldilut(float, float, float, float, float, float, float, float);
		float fphis(float phi, int dir=0);
		float fpsis(float phi, int dir=0);
		float fthetas(float theta, int dir=0);
		int spin();
		float poltrans_elg (float, float, int);//, int);
		float fepsilon(float, float, float);
	private:

};

#endif



