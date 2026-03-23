#define RadProcess_cxx
#include "RadProcess.h"
#include <limits>
#include <unistd.h>
using namespace std;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

double RadProcess::thetag_brem(double Ei, double Ef, double ZZ){
        Ei *= 1./m_el; Ef*= 1./m_el; // cf egs5 p 56 and paper from Bethe 1958
        double R = Ef/Ei;
        double Xi = (rand() /(double)RAND_MAX);
        double Th = Thsample(Ei);
        double Y = Th*Ei;
        double X=Y*Y;
        double MX = pow( (1-R) / (2.*Ei*R) ,2) + pow((pow(ZZ,1./3.)/(111.*(X+1.))),2) ;
        double Mun = pow((1-R) / (2.*Ei*R) ,2) + pow((pow(ZZ,1./3.)/222.),2) ;
        double Mzero = pow((1-R) / (2.*Ei*R) ,2) + pow((pow(ZZ,1./3.)/111.),2) ;
        double GX = 2.*R-3.*(1+R*R)-(4.+log(MX)) * ((1+R*R) - 4*X*R/pow(X+1,2));
        double Gzero = 2.*R-3.*(1+R*R)-(4.+log(Mzero)) * (1+R*R) ;
        double Gun = 2.*R-3.*(1+R*R)-(4.+log(Mun))  * (1+R*R - R);

        double NRinv = maximum(GX,Gzero,Gun);

        double Gtest = GX / NRinv;
        if (Xi<=Gtest) return Th;
        else return 0;

}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::Eg_rdint(float Eel, float QQ){
        double mel=0.00051099891;
        float R=(rand() /(double)RAND_MAX);
        float d= (alphaEM/PI *(log(QQ/pow(mel,2))-1 ));
        float Eg= Eel*pow(R,2./d);
        if (Eg>0.05 || (Eg/Eel)<0.0001) return 0; // exact cutoff IR: p 38 egs5
        else return Eg;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::EqRad_lenght(float QQ){
	double mel=0.00051099891;
	float d= 3./4.*(alphaEM/PI *(log(QQ/pow(mel,2))-1 ));
	return d;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::Eg_brem(float Eel,float Ecut,float d){
	float R=0;
	float Eg=0;
	for (int i=0;i<10000;i++){
		R=(rand() /(double)RAND_MAX)*0.3+0.7;
		Eg = Eel*pow(R,(1./(d*4./3.)));
		if (Eg>Eel*0.0001 && Eg<Ecut) return Eg;
	}
	//TF1 *f1=new TF1("f1",Form("%f*(1./x)*( (4./3.-4./3.*(x/%f)-pow((x/%f),2)) )",d,Eel,Eel), Eel*0.0001, Ecut);
	//float Eg= f1->GetRandom();
	return 0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

TLorentzVector RadProcess::TLBrDeviate(TLorentzVector el_in, float Z, float Ecut, float LL){
	// correct function: add angular rotation back. approx that theta=0 here valid for soft photons.
	float ein=el_in.E();
	float Eg=0, theta_gamma=0, phi1=0;
	for (int i=0;i<10;i++){
		Eg=Eg_brem(ein,Ecut,LL);
		if ((ein-Eg)<1.5) {
			continue;
		} else if (Eg>2){
			 // angular correction only for hard photons
			phi1=(rand() /(double)RAND_MAX)*2*PI;
			theta_gamma = thetag_brem(ein, ein-Eg,Z);
			el_in.SetPxPyPzE(el_in.Px()-Eg*sin(theta_gamma)*cos(phi1),el_in.Py()-Eg*sin(theta_gamma)*sin(phi1),el_in.Pz()-Eg*cos(theta_gamma),ein-Eg); 
		} else {
			el_in.SetPxPyPzE(el_in.Px(),el_in.Py(),el_in.Pz()-Eg,ein-Eg);
			break;
		}
	}
	//float phi1=(rand() /(double)RAND_MAX)*2*PI;
        //float theta_gamma = thetag_brem(ein, ein-Eg,Z); // very small for soft photon
	//el_in.SetPxPyPzE(el_in.Px()-Eg*sin(theta_gamma)*cos(phi1),el_in.Py()-Eg*sin(theta_gamma)*sin(phi1),el_in.Pz()-Eg*cos(theta_gamma),ein-Eg); 	
	//el_in.SetPxPyPzE(el_in.Px(),el_in.Py(),el_in.Pz()-Eg,ein-Eg); 	
	return el_in;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

double RadProcess::Thsample(double Einorm){

        double Xi = (rand() /(double)RAND_MAX);
        double Th = 1./Einorm*sqrt( Xi/ (1.-Xi+ 1./(PI*PI*Einorm*Einorm) ) );
        return Th;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

double RadProcess::VFluxFactor(double fQ2_max, double E, double nu){

        double fy = nu/E ;
        double me =510.998910e-6;
        double fQ2_min = me*me*fy*fy/(1 - fy);
        double result= (1/E) * (1./(137.036)) /(PI*fy)*( (1 - fy + fy*fy/2)*log(fQ2_max/fQ2_min) - (1 - fy));

        return result;

}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

float RadProcess::UNX(double AA, double ZZ){

	float Lrad, Lprad, fz, prefact, aa, unx0, rhoH=1;

        if (AA<1000){

                if (ZZ==1.){ // H or D
                        Lrad=5.31;
                        Lprad=6.144;
                        if (AA==1) rhoH=70.9*0.001;
			if (AA==2) rhoH=0.169;
                } else if (ZZ==2.){ // He
                        Lrad=4.79;
                        Lprad=5.621;
                        rhoH = 0.125;
                } else if (ZZ==3.) {
                        Lrad=4.74;
                        Lprad=5.805;
                        rhoH = 0.534 ; // Li
                } else if (ZZ==4.) {
                        Lrad=4.71;
                        Lprad=5.924;
                        rhoH = 1.848;
                } else
                {
                        Lrad=log(184.15*pow(ZZ,(-1./3.)));
                        Lprad=log(1194.*pow(ZZ,(-2./3.)));
                        if (ZZ==7) rhoH=0.682; // NH3 
                        else if (ZZ==79) rhoH = 1.932; // Au
                }
                prefact=1./716.408*AA;
                aa=alphaEM*ZZ;

                fz=aa*aa* ( (1./ (1.+aa*aa)) + 0.20206 - 0.0369*aa*aa+0.0083*pow(aa,4)-0.002*pow(aa,6));
                unx0=prefact*( ZZ*ZZ*(Lrad-fz)+ZZ*Lprad )*rhoH ;

        } else {
		unx0= 1./RadLenght[(int) (AA-1000)]; 
        }
	
	return unx0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

double RadProcess::BremstrahlungSpectraNgammadiff(double E,double nu,double ZZ,double AA, double dtarget){
        double fyy=nu/E;
	float unx0=UNX(AA,ZZ);
        return dtarget/2.*unx0*(1./nu)*(4./3. -4.*fyy/3. + pow(fyy,2));
	// flux in dN/dE GeV-1
        // *deltaE to get N gamma

}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::InitBrmProfile(float Eel, float Emin, float Emax){
	// function to normalize a flat generated distribution in Egamma to an actual bremstrahlung spectrum
	// called for first event, return the integral of N generated bremsstrahlung photons within Emin, Emax
	TF1 *f1 = new TF1 ("f1",Form("(1./x*(4./3.- (4.*x)/(3.*%f) + pow(x/%f,2) ))",Eel,Eel),Emin, Emax);
	float iB = f1->Integral(Emin, Emax);
	return iB;	
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::IntegralNph(float E, float kmin, float kmax, float d){

	float Ngamma;
	if (kmin<E/10000.){
                kmin=E/10000.;
        }
	// d in unit of rad lenght. (d*0.5/X0)
	Ngamma = d*(4./3.*log(kmax/kmin)-4./(3.*E)*(kmax-kmin)+(pow(kmax,2)-pow(kmin,2))/(2.*pow(E,2) ));
	return Ngamma;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::SampleNph(float Nph_steps[],float E, float AA, float ZZ, float d){
	// divide the bremsstrahlung full spectra into 10 steps and calculate the number of photons / electron within rad lenght (half target d)
	// adjusted to have about same number of photons/step. tested LH2 10cm
	// initialize only once
	float kmin, kmax;
	d=d*0.5;
	for (int i=0;i<10;i++){
		kmax = pow(2.5,i+1)*E/10000.; 
		kmin = pow(2.5,i)*E/10000.;
		Nph_steps[i]=IntegralNph(E, kmin, kmax,  d);
	}
	return IntegralNph(E, E/10000., pow(2.5,10)*E/10000.,  d);	

}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

float RadProcess::IntegralBr(float Eel, float Emin, float Emax,float AA, float ZZ, float LL){
	// need to be adapted. use gaus better from brm formula // only Z=1 A=1 now
	
	float unx0 = UNX(AA,ZZ);
	TF1 *f1 = new TF1 ("f1",Form("(%f*%f*0.5/x*(4./3.- (4.*x)/(3.*%f) + pow(x/%f,2) ))", LL, unx0,Eel,Eel),Emin, Emax);
	float iB = f1->Integral(Emin, Emax);
	return iB;	
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
float RadProcess::BeamProfileRescale_bmr(float Egam, float Eel, float Emin, float Emax, float ib){
	// proba density function normalize to 1*Delta_Eg between A and B. Same integral as of generation for photoproduction. 
	// give Delta=1 for bremsstrahlung? 
	// It returns P(Eg) Delta_Eg = dsigma(brem)/dE / integral of brem function from A to B * Delta_Eg
	// Egam has to be between A and B (cutoff RC or photon beam energy limits)
	// will allocate larger weight for lower energy photons within A and B limits P(A)=max, P(B)=min
        float res=0;
	//float iB = 0.353976;//0.63253;
        // res= (iB/0.63253)*(Emax-Emin)
	res =(Emax-Emin)*1./ib*(1./Egam*(4./3.- (4.*Egam)/(3.*Eel) + pow(Egam/Eel,2) ));
	return res;
}



//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

TLorentzVector RadProcess::ElectronRC(TLorentzVector el_in, float Ecut, float ZZ, float LL, float intcut, int &rad_event) {

	float test=1;
	for (int i=0;i<10;i++){
		test=(rand() /(double)RAND_MAX);
		if (test>intcut) break; // no radiation for this event
		else {
			el_in=TLBrDeviate(el_in, ZZ, Ecut, LL);
			rad_event+=1;
			intcut-=1;
			if (intcut<0) break; else continue;
		}	
	}
	return el_in;

}


