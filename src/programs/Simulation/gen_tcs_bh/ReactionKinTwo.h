#ifndef ReactionKinTwo_h
#define ReactionKinTwo_h
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
#include "PartUtils.h"

class ReactionKinTwo {

   public:
	void rotate_outphoto_backlab(TLorentzVector &, double);
	void rotate_outelectro_backlab(TLorentzVector &, double, double, double);
	int get_gammaPout(TLorentzVector &, TLorentzVector &, double&, double&, double&, double &, double&, double&, double, double, double, double, double, double , double);
	int get_gammaNout(TLorentzVector &, TLorentzVector &, double&, double&, double&, double &, double&, double&, double, double, double, double, double, double , double);
   private:

};

#endif



