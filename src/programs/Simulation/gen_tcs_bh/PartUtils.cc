#define PartUtils_cxx
#include "PartUtils.h"
using namespace std;

//-----------------------------------------------------------------------------------------------//

void PartUtils::Boost(TLorentzVector &LV, float betta, float gamma){
	LV.SetPxPyPzE((float)  LV.Px(),(float) LV.Py(),(float) gamma*(LV.Pz()-betta* LV.E()), (float) gamma*(LV.E()-betta*LV.Pz()) );
	return;
} 

void PartUtils::BoostBack(TLorentzVector &LV, float betta, float gamma){
	LV.SetPxPyPzE((float)  LV.Px(),(float) LV.Py(),(float) gamma*(LV.Pz()+betta* LV.E()), (float) gamma*(LV.E()+betta*LV.Pz()) );
	return;
} 

TLorentzVector PartUtils::LVBoost(TLorentzVector LV, float betta, float gamma){
	LV.SetPxPyPzE((float)  LV.Px(),(float) LV.Py(),(float) gamma*(LV.Pz()-betta* LV.E()), (float) gamma*(LV.E()-betta*LV.Pz()) );
	return LV;
} 

TLorentzVector PartUtils::LVBoostBack(TLorentzVector LV, float betta, float gamma){
	LV.SetPxPyPzE((float)  LV.Px(),(float) LV.Py(),(float) gamma*(LV.Pz()+betta* LV.E()), (float) gamma*(LV.E()+betta*LV.Pz()) );
	return LV;
}

TLorentzVector PartUtils::LVRot_direct_Z(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)-LV.Py()*sin(phi), LV.Px()*sin(phi)+ LV.Py()*cos(phi), LV.Pz(), LV.E());	
	return LV;
}
 
TLorentzVector PartUtils::LVRot_clock_Z(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)+LV.Py()*sin(phi), -LV.Px()*sin(phi)+ LV.Py()*cos(phi), LV.Pz(), LV.E());	
	return LV;
}

TLorentzVector PartUtils::LVRot_direct_Y(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)-LV.Pz()*sin(phi), LV.Py(), LV.Px()*sin(phi)+ LV.Pz()*cos(phi), LV.E());	
	return LV;
}

TLorentzVector PartUtils::LVRot_clock_Y(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)+LV.Pz()*sin(phi), LV.Py(), LV.Px()*sin(phi)- LV.Pz()*cos(phi), LV.E());	
	return LV;
}
 
TLorentzVector PartUtils::LVRot_direct_X(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px(), LV.Py()*cos(phi)-LV.Pz()*sin(phi), LV.Py()*sin(phi)+ LV.Pz()*cos(phi), LV.E());	
	return LV;
}

TLorentzVector PartUtils::LVRot_clock_X(TLorentzVector LV, float phi){
	LV.SetPxPyPzE(LV.Px(), LV.Py()*cos(phi)+LV.Pz()*sin(phi), LV.Py()*sin(phi)- LV.Pz()*cos(phi), LV.E());	
	return LV;
}

void PartUtils::Rot_direct_Z(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)-LV.Py()*sin(phi), LV.Px()*sin(phi)+ LV.Py()*cos(phi), LV.Pz(), LV.E());	
	return;
}
 
void PartUtils::Rot_clock_Z(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)+LV.Py()*sin(phi), -LV.Px()*sin(phi)+ LV.Py()*cos(phi), LV.Pz(), LV.E());	
	return;
}

void PartUtils::Rot_direct_Y(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)-LV.Pz()*sin(phi), LV.Py(), LV.Px()*sin(phi)+ LV.Pz()*cos(phi), LV.E());	
	return;
}

void PartUtils::Rot_clock_Y(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px()*cos(phi)+LV.Pz()*sin(phi), LV.Py(), LV.Px()*sin(phi)- LV.Pz()*cos(phi), LV.E());	
	return;
}
 
void PartUtils::Rot_direct_X(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px(), LV.Py()*cos(phi)-LV.Pz()*sin(phi), LV.Py()*sin(phi)+ LV.Pz()*cos(phi), LV.E());	
	return;
}

void PartUtils::Rot_clock_X(TLorentzVector &LV, float phi){
	LV.SetPxPyPzE(LV.Px(), LV.Py()*cos(phi)+LV.Pz()*sin(phi), LV.Py()*sin(phi)- LV.Pz()*cos(phi), LV.E());	
	return;
}

void PartUtils::FillArray_LV(TLorentzVector LV, double *Ar){
//void PartUtils::FillArray_LV(TLorentzVector LV, double (&Ar)[4]){
	Ar[0]=LV.E(); Ar[1]=LV.Px(); Ar[2]=LV.Py(); Ar[3]=LV.Pz();
	return;
} 

float PT(float px, float py){
	return sqrt(pow(px,2)+pow(py,2));
}

//-----------------------------------------------------------------------------------------------//
