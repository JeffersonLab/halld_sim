#ifndef CUDA_DBLREGGEHELPER_FASTETA
#define CUDA_DBLREGGEHELPER_FASTETA

#include "GPUManager/GPUCustomTypes.h"
#include <stdio.h>

static __device__ WCUComplex CHGM(double A, double B, double X) {
	double A0=A, A1=A, X0=X, HG = 0.0;
	double TBA, TB, TA, Y0=0.0, Y1=0.0, RG, LA = (int) A, NL, R, M, INF = pow(10.0,300.0);
	double sum1, sum2, R1, R2, HG1, HG2;
	if (B == 0.0 || B == -abs( (int) B)){
		HG = INF;
	} else if(A == 0.0 || X == 0.0) {
		HG = 1.0;
	} else if(A == -1.0){
		HG = 1.0 - X/B;
	} else if(A == B){
		HG = exp(X);
	} else if (A-B == 1.0){
		HG = (1.0+X/B)*exp(X);
	} else if (A == 1.0 && B == 2.0){
		HG = (exp(X)-1.0)/X;
	} else if(A == (int)A && A < 0.0){
		M = (int) -A;
		R = 1.0;
		HG = 1.0;
		for (int k = 1; k<= M ; k++) {
			R = R*(A+k-1.0)/k/(B+k-1.0)*X;
			HG+=R;
		}
	}
	if(HG != 0){
		WCUComplex CHG = {HG,0};
		return CHG;
	}

	if(X<0.0){
		A = B-A;
		A0 = A;
		X = fabs(X);
	}
	if(A<2.0) {NL = 0;}
	else{
		NL = 1;
		LA = (int) A;
		A  = A-LA-1.0;
	}
	for (int n = 0; n<= NL; n++) {
		if(A0 >= 2.0 ) { A+=1.0; }
		if(X <= 30.0 + fabs(B) || A < 0.0){
			HG = 1.0;
			RG = 1.0;
			for (int j = 1; j<= 500; j++) {
				RG = RG*(A+j-1)/(j*(B+j-1))*X;
				HG += RG;
				if(fabs(RG/HG) < pow(10.,-15.)) {
					if(n==0) {Y0 = HG;}
					if(n==1) {Y1 = HG;}
				}
				continue;
			}
		} else {
			TA = tgamma(A);
			TB = tgamma(B);
			TBA = tgamma(B-A);
			sum1 = 1.0;
			sum2 = 1.0;
			R1 = 1.0;
			R2 = 1.0;
			for (int i = 1; i<=8; i++) {
				R1 = - R1*(A+i-1)*(A-B+i)/(X*i);
				R2 = - R2*(B-A+i-1)*(A-i)/(X*i);
				sum1+=R1;
				sum2+=R2;
			}
			HG1 = TB/TBA*pow(X,-A)*cos(M_PI*A)*sum1;
			HG2 = TB/TA*exp(X)*pow(X,A-B)*sum2;
			HG = HG1+HG2;
		}
		if(n==0) {Y0 = HG;}
		if(n==1) {Y1 = HG;}
	}
	if(A0 >= 2.0){
		for (int i=1; i<=LA-1; i++) {
			HG = ((2.*A-B+X)*Y1+(B-A)*Y0)/A;
			Y0 = Y1;
			Y1 = HG;
			A += 1.;
		}
	}
	if(X0<0.0) {HG = HG*exp(X0);}
	A = A1;
	X = X0;

	WCUComplex CHG = {HG, 0};
	return CHG;
}


static __device__ WCUComplex cgamma(WCUComplex z,int OPT) {
	WCUComplex ui= {0,1};
	WCUComplex g, infini= 1e308+ 0.0*ui; // z0,z1
	double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
	double na=0.0,t,x1 = 1,y1=0.0,sr,si;
	int j,k;


	double a[] = {
		8.333333333333333e-02,
		-2.777777777777778e-03,
		7.936507936507937e-04,
		-5.952380952380952e-04,
		8.417508417508418e-04,
		-1.917526917526918e-03,
		6.410256410256410e-03,
		-2.955065359477124e-02,
		1.796443723688307e-01,
		-1.39243221690590};

	x = z.Re();
	y = z.Im();

	if (x > 171) return infini;
	if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
		return infini;
	else if (x < 0.0) {
		x1 = x;
		y1 = y;
		x = -x;
		y = -y;
	}
	x0 = x;
	if (x <= 7.0) {
		na = (int)(7.0-x);
		x0 = x+na;
	}
	q1 = sqrt(x0*x0+y*y);
	th = atan(y/x0);
	gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
	gi = th*(x0-0.5)+y*log(q1)-y;
	for (k=0;k<10;k++){
		t = pow(q1,-1.0-2.0*k);
		gr += (a[k]*t*cos((2.0*k+1.0)*th));
		gi -= (a[k]*t*sin((2.0*k+1.0)*th));
	}
	if (x <= 7.0) {
		gr1 = 0.0;
		gi1 = 0.0;
		for (j=0;j<na;j++) {
			gr1 += (0.5*log((x+j)*(x+j)+y*y));
			gi1 += atan(y/(x+j));
		}
		gr -= gr1;
		gi -= gi1;
	}


	if (x1 <= 0.0) {
		q1 = sqrt(x*x+y*y);
		th1 = atan(y/x);
		sr = -sin(M_PI*x)*cosh(M_PI*y);
		si = -cos(M_PI*x)*sinh(M_PI*y);
		q2 = sqrt(sr*sr+si*si);
		th2 = atan(si/sr);
		if (sr < 0.0) th2 += M_PI;
		gr = log(M_PI/(q1*q2))-gr;
		gi = -th1-th2-gi;
		x = x1;
		y = y1;
	}

	if (OPT == 0) {
		g0 = exp(gr);
		gr = g0*cos(gi);
		gi = g0*sin(gi);
	}
	g = gr + ui*gi;

	return g;
}



static __device__ WCUComplex V12(double alp1, double alp2, double eta) {

	WCUComplex zero = {0,0};
	if(alp1==alp2 ){return zero;}
	WCUComplex res = CHGM(-alp1, 1.-alp1+alp2, -1/eta);

	WCUComplex arg1 = {alp1-alp2,0};
	WCUComplex arg2 = {-alp2, 0};
	res *= cgamma(arg1,0)/cgamma(arg2,0);

	return res;
}



static __device__ WCUComplex GPU_DoubleRegge(int tau[2], GDouble s, GDouble si[2], GDouble alp[2], GDouble S0){

	WCUComplex ui=  {0,1};
	// signature factors:

	GDouble real,realE,imagE, imag;
	real = (-1.*ui*M_PI*alp[0]).Re();
	imag = (-1.*ui*M_PI*alp[0]).Im();

	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step = {realE, imagE};
	WCUComplex x0  = 1/2.*((double)tau[0] +step);

	real = (-1.*ui*M_PI*alp[1]).Re();
	imag = (-1.*ui*M_PI*alp[1]).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step1 = {realE, imagE};

	WCUComplex x1  = 1/2.*((double)tau[1] + step1);

	real = (-1.*ui*M_PI*(alp[0] - alp[1])).Re();
	imag = (-1.*ui*M_PI*(alp[0]-alp[1])).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step2 = {realE, imagE};

	WCUComplex x01 = 1/2.*((double)tau[0]*tau[1] + step2);

	real = (-1.*ui*M_PI*(alp[1] - alp[0])).Re();
	imag = (-1.*ui*M_PI*(alp[1]-alp[0])).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step3 = {realE, imagE};
	WCUComplex x10 = 1/2.*((double)tau[1]*tau[0] + step3);
	// double Regge vertices:

	double eta = S0*s/(si[0]*si[1]);
	WCUComplex V0 = V12(alp[0], alp[1], eta);
	WCUComplex V1 = V12(alp[1], alp[0], eta);

	GDouble up1 = pow(s/S0,alp[1])*pow(si[0]/S0,alp[0]-alp[1]);
	GDouble up2 = pow(s/S0,alp[0])*pow(si[1]/S0,alp[1]-alp[0]);

	// combine pieces:


	WCUComplex  t1 =up1*x1*x01*V1;
	WCUComplex t0 = up2*x0*x10*V0;

	WCUComplex arg0 = {-alp[0], 0};
	WCUComplex arg1 = {-alp[1], 0};
	return (t0+t1)*cgamma(arg0,0)*cgamma(arg1,0);;


}


static __device__ WCUComplex GPU_ampEtaPi0(GDouble par, int hel[3], GDouble inv[5], GDouble mass2[4], GDouble b_eta, GDouble S0, int charge ){

	WCUComplex zero =  {0,0};

        if(abs(hel[0]*hel[1]*hel[2]) != 1 || hel[0]==0){return zero;}

	GDouble s, s12,s23,t1,u3;
	GDouble m12, m22, m32, ma2;
	s   = inv[0];   s12 = inv[1];   s23 = inv[2];   t1  = inv[3];   u3  = inv[4];
	ma2 = mass2[0]; m12 = mass2[1]; m22 = mass2[2]; m32 = mass2[3];
	GDouble t2  = -t1+u3-s12+ma2+m12+m22;
	GDouble s13 = s-s12-s23+m12+m22+m32;

	// scalar part
	GDouble app = 0.9;     // slope of Regge trajectories alpha'
	GDouble alp0eta = app*t1 + 0.5;
	GDouble alp0pi0 = app*t2 + 0.5;
	GDouble alp1    = app*u3 + 0.5;

//	int tau[2] = {-1, -1};    // only vector exchange
	        int tau[2];

if(charge ==0){
        tau[0] = -1;
        tau[1] = -1;    // only vector exchange
}
else if(charge == 1){
        tau[0] = 1;
        tau[1]= -1; //for charged channel, a2 exchange?
}
	GDouble si[2]; 
	GDouble alp[2];
	si[0] = s12; si[1] = s23;
	alp[0] = alp0eta; alp[1] = alp1;

	WCUComplex ADR1 = GPU_DoubleRegge(tau, s, si, alp,S0); // fast eta
	GDouble fac1 =  sqrt(-t1/mass2[2]);
	GDouble fac3 = pow(-u3/4./mass2[0],abs((hel[1]-hel[2])/4.)); // hel[1,2] are twice the nucleon helicities!
	GDouble parity = pow(-1.0,(hel[1]-hel[2])/2.);
	if(hel[1] == -1){fac3 = fac3*parity;}

	GDouble Bot1 = exp(abs(b_eta)*t1);

	WCUComplex finalFactor = fac3*(Bot1*fac1*ADR1);
	return finalFactor;
}




static __device__ WCUComplex GPU_calcAmplitude(GDouble s, GDouble s12, GDouble s23, GDouble t1, GDouble u3, GDouble S0, GDouble b_eta, GDouble beamM2, GDouble p1M2, GDouble p2M2, GDouble recoilM2, int charge){

	GDouble param = b_eta;
	GDouble inv[5] = {s, s12, s23, t1 ,u3};
	GDouble mass2[4] = {beamM2, p1M2, p2M2, recoilM2};

	int hel[3] = {1,-1,-1};

	WCUComplex amp = GPU_ampEtaPi0(param, hel, inv, mass2, b_eta, S0, charge);


	return amp;
}




#endif
