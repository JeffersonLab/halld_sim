#define ReactionKinTwo_cxx
#include "ReactionKinTwo.h"
#include <limits>
#include <unistd.h>
using namespace std;

int ReactionKinTwo::get_gammaPout(
			TLorentzVector &LV_gout, TLorentzVector &LV_pout, double &Beta, double &Gamma, double &WW, double &costhcm, 
			double &EVirtual_CMeP, double &PVirtual_CMeP,  
			double Eb, double Eg_in, double Q2_in, double Q2_out, double TT, 
			double thmin=0, double thmax=180.
			){
	// this function calculates kinematic in CM frame and return virtual/real photon out + recoil nucleon
	// return 0 if kinematic not allowed or not within kin cuts
	double Eg_out, Pg_out, costhlab, PinCM, PoutCM, Thetagg_CMeP; 
	Kinematics kin; 
	PartUtils part_op;
	
	// lab kin, invariants
	WW=  kin.Wsqr(Eb, Eg_in, Q2_in);
	Eg_out = kin.fEg_out_lab(WW,TT,Q2_in);
        if (Eg_out<sqrt(Q2_out)) return 0;
        Pg_out = kin.fmom(Eg_out, sqrt(Q2_out));
        //tm = kin.ftmin(Eg_in,Eg_out ,-Q2_in,Q2_out);                 
        //if (TT>tm) return 0;
        costhlab= kin.fCosTh(TT, Eg_in, Eg_out, sqrt(pow(Egamma,2)+Q2_in), Pg_out, -Q2_in, Q2_out);
	if (fabs(costhlab)>1) return 0;

	// CM frame
	Beta=sqrt(pow(Eg_in,2)+Q2_in)/(Eg_in+M_Nucleon);
        Gamma=(Eg_in+M_Nucleon)/sqrt(WW);
        PinCM=M_Nucleon*sqrt((pow(Eg_in,2)+Q2_in)/WW);
        PoutCM=sqrt( (pow((WW-Q2_out-M_Nucleon*M_Nucleon),2)-4.*Q2_out*pow(M_Nucleon,2))/(4.*WW));
        EVirtual_CMeP=(WW-pow(M_Nucleon,2)+Q2_out)/(2.*sqrt(WW)) ;
       	PVirtual_CMeP=PoutCM;// sqrt( 1/(4.*WW)*( pow(WW-Q2_out-pow(M_Nucleon,2),2)-4*pow(M_Nucleon,2)*Q2_out));
        costhcm = kin.fCosTh(TT, sqrt(pow(PinCM,2)-Q2_in),EVirtual_CMeP, PinCM, PoutCM ,-Q2_in,Q2_out);
	if (fabs(costhcm)>1) return 0;
        Thetagg_CMeP=acos(costhcm);
	//if ((reaction==3 || reaction==4 || reaction==13) && 
	if ((Thetagg_CMeP*180./PI<thmin || Thetagg_CMeP*180./PI>thmax)) return 0;

	// 4-vectors CM	
	LV_gout.SetPxPyPzE((float) PoutCM*sin(Thetagg_CMeP),(float) 0. ,(float) PoutCM*cos(Thetagg_CMeP),sqrt(Q2_out+pow(PoutCM,2)));
        LV_pout.SetPxPyPzE((float) -PoutCM*sin(Thetagg_CMeP),(float) 0.,(float) -PoutCM*cos(Thetagg_CMeP),sqrt(pow(M_Nucleon,2)+pow(PoutCM,2)));

	// 4-vectors back to lab
	part_op.BoostBack(LV_gout, Beta, Gamma);
        part_op.BoostBack(LV_pout, Beta, Gamma);
	
	return 1;
}

//---------------------------------------------------------------------------------------------------------------

int ReactionKinTwo::get_gammaNout(
			TLorentzVector &LV_gout, TLorentzVector &LV_pout, double &Beta, double &Gamma, double &WW, double &costhcm, 
			double &EVirtual_CMeP, double &PVirtual_CMeP,  
			double Eb, double Eg_in, double Q2_in, double Q2_out, double TT, 
			double thmin=0, double thmax=180.
			){
	// this function calculates kinematic in CM frame and return virtual/real photon out + recoil nucleon
	// return 0 if kinematic not allowed or not within kin cuts
	double Eg_out, Pg_out, costhlab, PinCM, PoutCM, Thetagg_CMeP; 
	Kinematics kin; 
	PartUtils part_op;
	
	// lab kin, invariants
	WW=  kin.Wsqr(Eb, Eg_in, Q2_in);
	Eg_out = kin.fEg_out_lab(WW,TT,Q2_in);
        if (Eg_out<sqrt(Q2_out)) return 0;
        Pg_out = kin.fmom(Eg_out, sqrt(Q2_out));
        //tm = kin.ftmin(Eg_in,Eg_out ,-Q2_in,Q2_out);                 
        //if (TT>tm) return 0;
        costhlab= kin.fCosTh(TT, Eg_in, Eg_out, sqrt(pow(Egamma,2)+Q2_in), Pg_out, -Q2_in, Q2_out);
	if (fabs(costhlab)>1) return 0;

	// CM frame
	Beta=sqrt(pow(Eg_in,2)+Q2_in)/(Eg_in+M_Proton);
        Gamma=(Eg_in+M_Proton)/sqrt(WW);
        PinCM=M_Proton*sqrt((pow(Eg_in,2)+Q2_in)/WW);
        PoutCM=sqrt( (pow((WW-Q2_out-M_Neutron*M_Neutron),2)-4.*Q2_out*pow(M_Neutron,2))/(4.*WW));
        EVirtual_CMeP=(WW-pow(M_Neutron,2)+Q2_out)/(2.*sqrt(WW)) ;
        PVirtual_CMeP= PoutCM;//sqrt( 1/(4.*WW)*( pow(WW-Q2_out-pow(M_Neutron,2),2)-4*pow(M_Neutron,2)*Q2_out));
        costhcm = kin.fCosTh(TT, sqrt(pow(PinCM,2)-Q2_in),EVirtual_CMeP, PinCM, PoutCM ,-Q2_in,Q2_out);
	if (fabs(costhcm)>1) return 0;
        Thetagg_CMeP=acos(costhcm);
	//if ((reaction==3 || reaction==4 || reaction==13) && 
	if ((Thetagg_CMeP*180./PI<thmin || Thetagg_CMeP*180./PI>thmax)) return 0;

	// 4-vectors CM	
	LV_gout.SetPxPyPzE((float) PoutCM*sin(Thetagg_CMeP),(float) 0. ,(float) PoutCM*cos(Thetagg_CMeP),sqrt(Q2_out+pow(PoutCM,2)));
        LV_pout.SetPxPyPzE((float) -PoutCM*sin(Thetagg_CMeP),(float) 0.,(float) -PoutCM*cos(Thetagg_CMeP),sqrt(pow(M_Neutron,2)+pow(PoutCM,2)));

	// 4-vectors back to lab
	part_op.BoostBack(LV_gout, Beta, Gamma);
        part_op.BoostBack(LV_pout, Beta, Gamma);
	
	return 1;
}
		
//---------------------------------------------------------------------------------	
//---------------------------------------------------------------------------------	

void ReactionKinTwo::rotate_outphoto_backlab(TLorentzVector &LV, double phib){
	PartUtils pu; 
	pu.Rot_clock_Z(LV,(phib));
	return;
}
//---------------------------------------------------------------------------------	
//---------------------------------------------------------------------------------	

void ReactionKinTwo::rotate_outelectro_backlab(TLorentzVector &LV, double phib, double thetab, double phip){
	PartUtils pu;
	pu.Rot_direct_Z(LV, phip);
	pu.Rot_direct_Y(LV,thetab);
	pu.Rot_direct_Z(LV,(phib+PI));
	return;
}

