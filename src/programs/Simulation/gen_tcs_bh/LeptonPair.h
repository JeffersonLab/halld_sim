#ifndef LeptonPair_h
#define LeptonPair_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "RadProcess.h"
#include "Utils.h"
#include <unistd.h>
#include "Constants.h"
#include "Kinematics.h"

class LeptonPair {

   public:
	void get_pair(TLorentzVector &LV_minus, TLorentzVector &LV_plus, float Epart, float Ppart, float Thpart, float thetaCM, float phiCM, float lepm);

   private:

};

#endif



