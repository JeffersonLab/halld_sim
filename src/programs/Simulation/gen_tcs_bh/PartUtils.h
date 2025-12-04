#ifndef PartUtils_h
#define PartUtils_h
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include "Constants.h"
#include "TROOT.h"
#include "TLorentzVector.h"

class PartUtils {

	public:
		void Boost(TLorentzVector &LV, float betta, float gamma);
		void BoostBack(TLorentzVector &LV, float betta, float gamma);
		TLorentzVector LVBoost(TLorentzVector LV, float betta, float gamma);
		TLorentzVector LVBoostBack(TLorentzVector LV, float betta, float gamma);
		float PT(float,float);
		void Rot_direct_X(TLorentzVector&,float);
		void Rot_direct_Y(TLorentzVector&,float);
		void Rot_direct_Z(TLorentzVector&,float);
		void Rot_clock_X(TLorentzVector&,float);
		void Rot_clock_Y(TLorentzVector&,float);
		void Rot_clock_Z(TLorentzVector&,float);
		TLorentzVector LVRot_direct_X(TLorentzVector,float);
		TLorentzVector LVRot_direct_Y(TLorentzVector,float);
		TLorentzVector LVRot_direct_Z(TLorentzVector,float);
		TLorentzVector LVRot_clock_X(TLorentzVector,float);
		TLorentzVector LVRot_clock_Y(TLorentzVector,float);
		TLorentzVector LVRot_clock_Z(TLorentzVector,float);
		void FillArray_LV(TLorentzVector LV, double *Ar);
		//void FillArray_LV(TLorentzVector LV, double (&Ar)[4]);
	private:
};

#endif


