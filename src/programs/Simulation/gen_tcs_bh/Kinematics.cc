#define Kinematics_cxx
#include <iostream>
#include <fstream>
#include "Kinematics.h"
#include <limits>
#include <unistd.h>
using namespace std;

float Kinematics::fQ2(float E, float Ep, float theta){
	return 4.*E*Ep*pow(sin(theta*0.5),2);
}

//---------------------------------------

float Kinematics::fXbj(float qq, float nu){
	return qq/(2.*M_Nucleon*nu);
}

//---------------------------------------

float Kinematics::fXbj(float E, float Ep, float theta){
	return fQ2(E, Ep, theta)/(2.*M_Nucleon*(E-Ep));
}

//---------------------------------------

float Kinematics::Wsqr(float E, float nu, float QQ){
	return pow(M_Nucleon,2)+2*M_Nucleon*nu -QQ;
}

//---------------------------------------------------

float Kinematics::fEg_out_lab(float Wsq, float t, float QQ){
	return (Wsq+t+QQ-pow(M_Nucleon,2))*0.5/M_Nucleon;
}

//-----------------------------------------------------

float Kinematics::fmom(float E, float M){
	return sqrt(pow(E,2)-pow(M,2));
}

//--------------------------------------------------------

float Kinematics::fen(float P, float M){
	return sqrt(pow(P,2)+pow(M,2));
}

//--------------------------------------------------------

float Kinematics::fmass(float E, float P){
	return sqrt(pow(E,2)-pow(P,2));
}

//--------------------------------------------------------

float Kinematics::ftmin(float Ein, float Eout, float M2in, float M2out){
	return (2.*fmom(Ein,sqrt(M2in))*fmom(Eout,sqrt(M2out)))-(2.*Ein*Eout-M2in-M2out);
}

//--------------------------------------------------------------------

float Kinematics::fCosTh(float t, float Ein, float Eout,float Pin, float Pout, float M2in, float M2out){
	return (t+2.*Ein*Eout-M2in-M2out)/(2.*Pin*Pout );
}





