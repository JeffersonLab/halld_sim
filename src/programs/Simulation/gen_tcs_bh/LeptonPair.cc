#define LeptonPair_cxx
#include "LeptonPair.h"
#include <limits>
#include <unistd.h>
using namespace std;


void LeptonPair::get_pair(TLorentzVector &LV_minus, TLorentzVector &LV_plus, float Epart, float Ppart, float Thpart, float thetaCM, float phiCM, float leptmass){

	float masspart, BetaV, GammaV, Plepton, Eplus_CM, Eminus_CM, Pzplus_CM, Pzminus_CM;
	//TLorentzVector LV_minus, LV_plus;

	// particle that decays into lepton pair. start in collision system CM
	masspart = sqrt(pow(Epart,2)-pow(Ppart,2));
	BetaV= Ppart/Epart;
        GammaV= Epart/masspart;
        Plepton=sqrt(pow(masspart,2)/4.-pow(leptmass,2));

	LV_minus.SetPxPyPzE(Plepton*sin(thetaCM)*cos(phiCM),Plepton*sin(thetaCM)*sin(phiCM),Plepton*cos(thetaCM),masspart/2.);
	LV_plus.SetPxPyPzE(-Plepton*sin(thetaCM)*cos(phiCM),-Plepton*sin(thetaCM)*sin(phiCM),-Plepton*cos(thetaCM),masspart/2.);

	// boost to CM eP (before rotation around y=y')
	Eplus_CM= GammaV*(sqrt(Qp2)/2.+BetaV*LV_plus.Pz());
        Eminus_CM= GammaV*(sqrt(Qp2)/2.+BetaV*LV_minus.Pz());
        Pzplus_CM= GammaV*(BetaV*sqrt(Qp2)/2.+LV_plus.Pz());
        Pzminus_CM= GammaV*(BetaV*sqrt(Qp2)/2.+LV_minus.Pz());
	
	// set and rotation
	LV_minus.SetPxPyPzE((float) Pzminus_CM*sin(Thpart)+LV_minus.Px()*cos(Thpart), LV_minus.Py(), (float) Pzminus_CM*cos(Thpart)-LV_minus.Px()*sin(Thpart), Eminus_CM);
        LV_plus.SetPxPyPzE((float) Pzplus_CM*sin(Thpart)+LV_plus.Px()*cos(Thpart), LV_plus.Py(), (float) Pzplus_CM*cos(Thpart)-LV_plus.Px()*sin(Thpart), Eplus_CM);	

	// back to lab (with gamma along z frame)
	//LV_minus.SetPxPyPzE((float)LV_minus.Px(),(float) LV_minus.Py(),(float) Gamma*(LV_minus.Pz()+Beta*LV_minus.E()), (float) Gamma*(LV_minus.E()+Beta*LV_minus.Pz()) );
        //LV_plus.SetPxPyPzE((float) LV_plus.Px(), (float) LV_plus.Py(), (float) Gamma*(LV_plus.Pz()+Beta*LV_plus.E()),   (float) Gamma*(LV_plus.E() +Beta*LV_plus.Pz()) );


	return;
}


