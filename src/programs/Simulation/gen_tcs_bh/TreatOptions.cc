#define TreatOptions_cxx
#include "TreatOptions.h"
#include <limits>
#include <unistd.h>
using namespace std;

/*
extern int NEB, NT ,NQp2, NTh, NPhi;
extern double EMINI, EMAXI, TMINI, TMAXI, QP2MINI, QP2MAXI, PHIMINI, PHIMAXI, THMINI, THMAXI; 
extern const double alphaEM, m_el, PI;

extern double yy_min, yy_max, Phi_LH, tt, Egamma, Qp2, Q2, Phi_CMV, Theta_CMV, yy, phi_beam, Xbj, Eproton, crossA, eta_min, eta_max;
extern int thmE, thmQ, thmT;
extern float thmEmin, thmEmax, thmQmin, thmQmax, thmTmin, thmTmax;
*/

 int NEB, NT ,NQp2, NTh, NPhi;
 int NPhis, NPsis; //added
 float CosThetaMaxTCS; // added
 double EMINI, EMAXI, TMINI, TMAXI, QP2MINI, QP2MAXI, PHIMINI, PHIMAXI, THMINI, THMAXI; 
 double PHISMINI, PHISMAXI, PSISMINI, PSISMAXI; //added
 extern const double alphaEM, m_el, PI;

 extern double yy_min, yy_max, Phi_LH, tt, Egamma, Qp2, Q2, Phi_CMV, Theta_CMV, yy, phi_beam, Xbj, Eproton, crossA, eta_min, eta_max;
 int thmE, thmQ, thmT;
 float thmEmin, thmEmax, thmQmin, thmQmax, thmTmin, thmTmax;

////////////////////////////////////////////////////


void settcsbins(int Qp2_LimitType3){
	if (Qp2_LimitType3==1){
	    // Binning for Low Q'2 TCS
	    NEB=13; NT=30; NQp2=15; NTh=25; NPhi=20; NPhis=20; NPsis=20; // tablev11 costh<1
	    CosThetaMaxTCS = 1;// 0.9995;
	    EMINI=5.0; EMAXI=11.5; TMINI=0.04; TMAXI=1.54; QP2MINI=0.8; QP2MAXI=5.3; // table v11
	    PHIMINI=3.; PHIMAXI=363.; THMINI=30.; THMAXI=155.; 
	    PHISMINI = 0; PHISMAXI=360; PSISMINI = 0; PSISMAXI = 360;

	    // table cuts (TCS)
	    thmE=14; thmQ=30; thmT=40;
	    thmEmin=5; thmEmax=12; thmQmin=0.8; thmQmax=9.3; thmTmin=0.02; thmTmax=2.02;
	    }

	// ****** High Q'2 TCS ******
	if (Qp2_LimitType3==2){
	    // Binning for High Q'2 TCS
	    NEB=13; NT=33; NQp2=18; NTh=20; NPhi=20; NPhis=20; NPsis=20; // tablev11 costh<1
	    CosThetaMaxTCS = 1;// 0.9995;
	    EMINI=5.0; EMAXI=11.5; TMINI=0.04; TMAXI=2.02; QP2MINI=3.8; QP2MAXI=9.2; // table v11
	    PHIMINI=3.; PHIMAXI=363.; THMINI=30.; THMAXI=150.; PHISMINI = 0; PHISMAXI=360; PSISMINI = 0; PSISMAXI = 360;

	    // table cuts (TCS)
	    thmE=14; thmQ=20; thmT=40;
	    thmEmin=5; thmEmax=12; thmQmin=3.8; thmQmax=9.3; thmTmin=0.02; thmTmax=2.02;
	    }

	// ****** Full Q'2 TCS ****** to be edited when table values becomes available
	if (Qp2_LimitType3==3){
	    // Binning for Full Q'2 TCS
	    NEB=13; NT=33; NQp2=18; NTh=20; NPhi=20; NPhis=20; NPsis=20; // tablev11 costh<1
	    CosThetaMaxTCS = 1;// 0.9995;
	    EMINI=5.0; EMAXI=11.5; TMINI=0.04; TMAXI=2.02; QP2MINI=3.8; QP2MAXI=9.2; // table v11
	    PHIMINI=3.; PHIMAXI=363.; THMINI=30.; THMAXI=150.; PHISMINI = 0; PHISMAXI=360; PSISMINI = 0; PSISMAXI = 360;

	    // table cuts (TCS)
	    thmE=14; thmQ=20; thmT=40;
	    thmEmin=5; thmEmax=12; thmQmin=0.8; thmQmax=9.3; thmTmin=0.02; thmTmax=2.02;
	    }
	return ;
}

///////////////////////////////////////////////////////////////////////////////////

int VariablesValues_PSTCSFIX(){
	if (Emax>Eelectron && beamtype>0) {
		Emax=Eelectron;
		cout<<"change Emax to "<<Eelectron<<endl;
	}
	if (Q2max>0.3 && beamtype>0) {
		cout<<"wrong value of Q2max: "<<Q2max<<endl;
		Q2max=0.3;
		cout<<"maximal value for quasi real photon: Q2=0.3 GeV2, changed"<<endl;
	}
	if (outlepton==1) {
		M_lepton= m_el;
		cout<<"Work with electron pair";
	}
	if (outlepton==2){
		M_lepton=M_Muon;
		cout<<"Work with muon pair";
	}
	M_Nucleon=M_Proton;
	if (beamtype==1){
		cout<<"Work with electron beam at "<<Eelectron<<"GeV"<<endl;
	}
	return 1;

}

///////////////////////////////////////////////////////////////////////////////////////////////

int VariablesValues_TCS(){

	cout<<"Table limits: "<<endl;
	cout<<"Egam: "<<EMINI<<" "<<EMAXI<<"   if el beam, E= "<<Eelectron<<endl;
	cout<<"-t: "<<TMINI<<" "<<TMAXI<<"   Q'2: "<<QP2MINI<<" "<<QP2MAXI<<endl;
	cout<<"ThCM: "<<THMINI<<" " <<THMAXI<<"   Phi: default= 0, 2pi\n"<<endl;

	if (Emax>Eelectron && beamtype>0) {
		Emax=Eelectron;
		cout<<"change Emax to "<<Eelectron<<endl;
	}
	if (Emin<EMINI){ 
		Emin=EMINI;
		cout<<"change Emin to "<<EMINI<<endl;
	}
	if (Emax>EMAXI){ 
		Emax=EMAXI; 
		cout<<"change Emax to "<<EMAXI<<endl;
	}
	if (mt_min<TMINI){ 
		mt_min=TMINI;
		cout<<"change -t min to "<<mt_min<<endl;
	} 
	if (mt_max>TMAXI){ 
		mt_max=TMAXI;
		cout<<"change -t max to "<<mt_max<<endl;
	}
	if (Qp2min<QP2MINI){ 
		cout<<"wrong value for Q'2 min: "<<Qp2min <<endl;
		Qp2min=QP2MINI;
		cout<<"change Q'2 min to "<<Qp2min<<endl;
	} 
	if (Qp2max>QP2MAXI){ 
		Qp2max=QP2MAXI;
		cout<<"change Q'2 max to "<<Qp2max<<endl;
	}
	if (Q2max>0.3 && beamtype>0) {
		cout<<"wrong value of Q2max: "<<Q2max<<endl;
		Q2max=0.3;
		cout<<"maximal value for quasi real photon: Q2=0.3 GeV2, changed"<<endl;
	}
	if (theta_min<THMINI) {
		theta_min=THMINI;
		cout<<"theta min changed to "<<THMINI<<endl;
	}

	if (theta_max>THMAXI){
		theta_max=THMAXI;
		cout<<"theta max changed to "<<THMAXI<<endl;
	}

	if (outlepton==1) {
		M_lepton= m_el;
		cout<<"Work with electron pair";
	}
	if (outlepton==2){
		M_lepton=M_Muon;
		cout<<"Work with muon pair";
	}
	if (protonorneutron==1) {
		cout<<" and with proton target"<<endl;
		M_Nucleon=M_Proton;
	} else if (protonorneutron==2) {
		M_Nucleon=M_Neutron;
		cout<<" and with neutron target"<<endl;
		cout<<"Warning: no table yet"<<endl;
	}
	if (beamtype==1){
		cout<<"Work with electron beam at "<<Eelectron<<"GeV"<<endl;
	}

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////

int RandomGen_TCS(){
	Egamma = (rand() /(double)RAND_MAX)* (Emax-Emin)+Emin ;
	Qp2=(rand() /(double)RAND_MAX)* (Qp2max-Qp2min)+ Qp2min ;
	tt =-((rand() /(double)RAND_MAX)* (mt_max-mt_min)+ mt_min );
	Phi_CMV= (rand() /(double)RAND_MAX)* (PHIMAXI-PHIMINI)*PI/180. + PHIMINI*PI/180.;
	//Phi_CMV= ((rand() /(double)RAND_MAX)* (PHIMAXI-PHIMINI) + PHIMINI)*PI/180.;
	Theta_CMV=  (rand() /(double)RAND_MAX)* ( theta_max-theta_min)*PI/180.+theta_min*PI/180.;
	//Theta_CMV=  ((rand() /(double)RAND_MAX)* ( theta_max-theta_min)+theta_min)*PI/180.;
	yy=Egamma/Eelectron;
	phi_beam= (rand() /(double)RAND_MAX)*2.*PI;
	if (beamtype==1) {
		Q2min = pow(510.998910e-6*yy,2)/(1 - yy);
		Q2=(rand() /(double)RAND_MAX)* Q2max ;
		if (Q2<Q2min) return 0;
		Phi_LH=(rand() / (double)RAND_MAX) *2.*PI;
	}
	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////

double L_BH_TCS(double eg, double q_out_th_lab, double pair_e_th_cm, double pair_e_phi_cm){

	double LL,a,b,c;
	double s,mt,ga_mom_out_lab,qin_cm,ega_in_cm, qout_cm, ega_out_cm, th_gaout_cm, r,deltaT,b1,b2,bbb,baba;
	double mel=0.00051099891;

	s = pow(M_Nucleon, 2.) + 2. * M_Nucleon * eg;

	a = pow(M_Nucleon + eg, 2.) - pow(eg * cos(q_out_th_lab), 2.);
	b = 0.5 * (s - pow(M_Nucleon,2.) + Qp2) * eg * cos(q_out_th_lab);
	c = pow(0.5 * (s - pow(M_Nucleon,2.) + Qp2), 2.) - Qp2 * pow(M_Nucleon + eg, 2.);

	if ((b*b + a*c) > 0){
		ga_mom_out_lab = 1./a * ( b + sqrt(pow(b,2.) + a * c));
	} else{
		printf("\n Error in kinematics outgoing photon lab momentum");
		return 0.;
	}

	mt = 2. * M_Nucleon * (eg - sqrt(Qp2 + pow(ga_mom_out_lab,2.)));
	qin_cm = sqrt((pow(s,2.)-2.*s*pow(M_Nucleon,2.)+pow(M_Nucleon,4.))/(4.*s));
	ega_in_cm = qin_cm;
	qout_cm = sqrt((pow(s,2.)-2.*s*(pow(M_Nucleon,2.)+Qp2)+ pow(Qp2-pow(M_Nucleon,2.),2.))/(4.*s));
	ega_out_cm = sqrt(pow(qout_cm, 2.) + Qp2);
	r = sqrt(pow(s-Qp2-pow(M_Nucleon,2.),2.)-4.*Qp2*pow(M_Nucleon,2.));
	th_gaout_cm = acos((-mt - Qp2 + 2. * ega_in_cm * ega_out_cm) / (2. * qin_cm * qout_cm));
	deltaT = sin(th_gaout_cm)/2./sqrt(s)*r;

	baba=sqrt(1.-4.*mel*mel/Qp2);
	b2 = 2.*(s-M_Nucleon*M_Nucleon)*sqrt(Qp2)*deltaT/r;
	b1 = sqrt(pow((Qp2-(-mt)),2) - pow(b2,2));
	bbb = baba*( b1*cos(pair_e_th_cm) - b2*sin(pair_e_th_cm)*cos(pair_e_phi_cm) );

	LL = (pow((Qp2+mt),2)-pow(bbb,2))/4.;
	return LL;

}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

double BH_anal_TCS_exact(double eg, double q_out_th_lab, double pair_e_th_cm, double pair_e_phi_cm){

	double result,calcul,LL,bet,AA,BB,b2,b1,bebe,aba;
	double s,  a, b, c, ga_mom_out_lab, mt, qin_cm,Egammaa_in_cm, qout_cm, Egammaa_out_cm, th_gaout_cm, r,  deltaT, f1, f2;
	double mel=0.00051099891;

	s = pow(M_Nucleon, 2.) + 2. * M_Nucleon * Egamma;
	a = pow(M_Nucleon + Egamma, 2.) - pow(Egamma * cos(q_out_th_lab), 2.);
	b = 0.5 * (s - pow(M_Nucleon,2.) + Qp2) * Egamma * cos(q_out_th_lab);
	c = pow(0.5 * (s - pow(M_Nucleon,2.) + Qp2), 2.) - Qp2 * pow(M_Nucleon + Egamma, 2.);

	if ((b*b + a*c) > 0){
		ga_mom_out_lab = 1./a * ( b + sqrt(pow(b,2.) + a * c));
	} else {
		printf("\n Error in kinematics outgoing photon lab momentum");
		return 0.;
	}

	mt = 2. * M_Nucleon * (Egamma - sqrt(Qp2 + pow(ga_mom_out_lab,2.)));

	qin_cm = sqrt((pow(s,2.)-2.*s*pow(M_Nucleon,2.)+pow(M_Nucleon,4.))/(4.*s));
	Egammaa_in_cm = qin_cm;
	qout_cm = sqrt((pow(s,2.)-2.*s*(pow(M_Nucleon,2.)+Qp2)+ pow(Qp2-pow(M_Nucleon,2.),2.))/(4.*s));
	Egammaa_out_cm = sqrt(pow(qout_cm, 2.) + Qp2);
  	th_gaout_cm = acos((-mt - Qp2 + 2. * Egammaa_in_cm * Egammaa_out_cm) / (2. * qin_cm * qout_cm));
	r=sqrt(pow(s-Qp2-pow(M_Nucleon,2.),2.)-4.*Qp2*pow(M_Nucleon,2.));
	deltaT=sin(th_gaout_cm)/2./sqrt(s)*r;

	f1 =FormFactors(1,1 ,mt);
		f2 =FormFactors(2,1,mt);

	LL = L_BH_TCS(Egamma,q_out_th_lab,pair_e_th_cm,pair_e_phi_cm);

	bet=sqrt(1.-4.*mel*mel/Qp2);
	b2 = 2.*(s-M_Nucleon*M_Nucleon)*sqrt(Qp2)*deltaT/r;
	b1 = sqrt(pow((Qp2+mt),2) - pow(b2,2));
	bebe= bet*( b1*cos(pair_e_th_cm) - b2*sin(pair_e_th_cm)*cos(pair_e_phi_cm) );
	aba=bet*r*cos(pair_e_th_cm);

	AA = pow((s - M_Nucleon*M_Nucleon),2)*pow(deltaT,2) + mt*aba*(aba+bebe) - pow((M_Nucleon*bebe),2)
		+ mt*(4.*M_Nucleon*M_Nucleon+mt)*Qp2
		+mel*mel/LL *( pow(( (Qp2+mt)* (aba+bebe) - (s - M_Nucleon*M_Nucleon)*bebe),2) - mt*(4.*M_Nucleon*M_Nucleon+mt)*pow((Qp2+mt),2));


	BB = pow((Qp2-mt),2) + pow(bebe,2) + 8.*pow(mel,2)*Qp2
		- 4.*pow(mel,2)*(-mt+2.*pow(mel,2))/LL*pow((Qp2+mt),2);

	calcul = (pow(f1,2)+mt/(4.*pow(M_Nucleon,2))*pow(f2,2))*AA/mt + pow((f1+f2),2)*BB/2.;
	result = pow(alphaEM*alphaEM/4./PI, 3.)*bet/(4.*PI*pow((s-pow(M_Nucleon,2)),2)*mt*LL)*calcul;
	result *= pow(.197, 2.) * pow(10., 7.);
	return result;

}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
