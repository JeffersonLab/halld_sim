#define Polar_cxx
#include <iostream>
#include <fstream>
#include "Polar.h"
#include "Constants.h"
#include <limits>
#include <unistd.h>
using namespace std;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::poltrans_elg(float polin, float y, int poltype){
	float pol=polin;

	if (poltype==0) pol= polin *y*(4.-y)/(4.-4.*y+3.*pow(y,2.)) ;

	return pol ; //polin *y*(4.-y)/(4.-4.*y+3.*pow(y,2.)) ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::fepsilon(float qq, float nu, float theta){
	return pow(1.+2*(pow(nu,2)+qq)/qq*pow(tan(theta/2.),2),-1);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int Polar::spin(){
	int dir;
	float test=(rand() /(double)RAND_MAX);
        if (test>=0.5) dir = 1; else dir = -1;
	return dir;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::fphis(float phi, int dir){

	float phis=0;
	if (dir==1) {
                phis= phi; 
	} else if (dir==2){
		phis=  phi -PI/2.;
	} 
	if (phis>=2*PI) phis -= 2*PI;
	if (phis<0) phis += 2.*PI;
	return phis;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::fthetas(float theta, int dir){

	if (dir==1 || dir==2){
		return PI/2.+theta;
	} else if (dir==3){
		return theta;
	} else return 0;

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::fpsis(float phi, int dir){

	float psis=0;
	if (dir==1) {
                psis= phi; 
	} else if (dir==2){
		psis=  phi -PI/2.;
	} else if (dir==3){
		psis=  phi -PI/4.;
	} 
	if (psis>=2*PI) psis -= 2*PI;
	if (psis<0) psis += 2.*PI;
	return psis;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::cross_poldilut(float p, float m, float spindir, float dil){
	float cross=(p+m)*0.5;
	if ( spindir>=0 ){
                        cross = p*dil + (1-dil)*cross;
        } else {
                        cross = m*dil + (1-dil)*cross;
        }
	return cross;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

float Polar::cross_doublepoldilut(float polpp, float polpm, float polmp, float polmm, float b_spindir, float t_spindir, float beamdil, float targetdil){
	float cross = (polpp+polpm+polmp+polmm)*0.25;
	if (t_spindir>=0) {
        	if (b_spindir>=0) {
                        cross = polpp*targetdil*beamdil+(1-targetdil*beamdil)*cross;
                } else {
                	cross = polmp*targetdil*beamdil+(1-targetdil*beamdil)*cross;
                }
        } else {
                if (b_spindir>=0) {
                        cross = polpm*targetdil*beamdil+(1-targetdil*beamdil)*cross;
                } else {
                        cross = polmm*targetdil*beamdil+(1-targetdil*beamdil)*cross;
                }
        }
	return cross;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------




