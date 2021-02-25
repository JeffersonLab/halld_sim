#include <cassert>
#include <complex>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "dblReggeMod.h"

dblReggeMod::dblReggeMod( const vector< string >& args ) :
	UserAmplitude< dblReggeMod >( args )
{ 
	assert( args.size() == 5 );
	//	j = atoi( args[0].c_str() );
	//	fast = atoi( args[1].c_str() );
	a_eta = AmpParameter(args[0]);
	a_pi = AmpParameter(args[1]);
	b_eta = AmpParameter( args[2] );
	b_pi = AmpParameter( args[3]);
	S0 = AmpParameter( args[4] );

	registerParameter( a_pi );
	registerParameter( a_eta );
	registerParameter( b_pi );
	registerParameter( b_eta );
	registerParameter( S0 );
}

complex< GDouble >
dblReggeMod::calcAmplitude( GDouble** pKin ) const {
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );

	TLorentzVector resonance = p1 + p2;

	double mass2[4] = {beam.M2(), p1.M2(), p2.M2(), recoil.M2()};

	double s12 = resonance.M2();
	//double s13 = (p1 + recoil).M2();
	double s23 =  (p2 + recoil).M2();
	double t1 = (beam - p1).M2();
	double t2 = (beam - p2).M2();
	double s = (recoil + p1 + p2).M2();
	double u3 = t1 + t2 +s12 - (beam.M2() + p1.M2() + p2.M2());


	double param[4] = {a_eta,a_pi,b_eta,b_pi};
	double inv[5] = {s, s12, s23, t1 ,u3};
//	double varGJ[5] = {s,s12,u3, theta, phi};

	//cout << "s: " << s << endl;
	//cout << "s12: " << s12 << endl;
	//cout << "s23: " << s23 << endl;
	//cout << "t1: " << t1 << endl;
	//cout << "u3: " << u3 << endl;

	int hel[3] = {1,-1,-1};
	std::complex<double> amp = ampEtaPi0(param, hel, inv, mass2);

//	if(abs(amp) > 1000)
//		cout << "Amplitude: " << abs(amp) << endl;

	return amp;
}


std::complex<double> dblReggeMod::ampEtaPi0(double par[4], int hel[3], double inv[5], double mass2[4]) const{

	std::complex<double> zero (0,0);
	// check that inputs are physical
	//double delta4 = GramDet4(mass2, inv);
	//if(delta4>0){return zero;}

	//	if(abs(hel[0]*hel[1]*hel[2]) != 1 || hel[0]==0){return zero;}


	double s,s12,s23,t1,u3;
	double m12, m22, m32, ma2, mb2;
	s   = inv[0];   s12 = inv[1];   s23 = inv[2];   t1  = inv[3];   u3  = inv[4];
	ma2 = mass2[0]; m12 = mass2[1]; m22 = mass2[2]; m32 = mass2[3];
	double t2  = -t1+u3-s12+ma2+m12+m22;
	double s13 = s-s12-s23+m12+m22+m32;

	// scalar part
	double app = 0.9;     // slope of Regge trajectories alpha'
	double alp0eta = app*t1 + 0.5;
	double alp0pi0 = app*t2 + 0.5;
	double alp1    = app*u3 + 0.5;

//	cout <<"t1, t2, u3: " << t1 <<", " <<t2<<", "<<u3 << endl;
//	cout <<"s, s12, s23: " << s <<", "<<s12<<", "<<s23 <<endl;


	int tau[2] = {-1, -1};    // only vector exchange
	double si[2] = {s12,s23};
	double alp[2] = {alp0eta, alp1};


	//cout << "fast: " << fast <<endl;
 	std::complex<double> ADR1 = DoubleRegge(tau, s, si, alp); // fast eta

	//cout <<"ADR1: " << ADR1 <<endl;

	si[1] = s13; alp[0] = alp0pi0;
	//	if(fast==2){
	std::complex<double> ADR2 = DoubleRegge(tau, s, si, alp); // fast pi0
	//	}
	// helicity part
	double fac1 =  sqrt(-t1/mass2[2]);
	double fac2 = sqrt(-t2/mass2[2]);    // use the pion mass in both fac1 and fac2
	double fac3 = pow(-u3/4./mass2[0],abs((hel[1]-hel[2])/4.)); // hel[1,2] are twice the nucleon helicities!
	double parity = pow(-1,(hel[1]-hel[2])/2.);
	if(hel[1] == -1){fac3 = fac3*parity;}

	double Bot1 = a_eta*exp(b_eta*t1);
	double Bot2 = a_pi*exp(b_pi*t2);
//	cout <<"fac3: " << fac3 <<endl;
//	cout <<"fac2: " << fac2 <<endl;
//	cout <<"fac1: " << fac1 <<endl;
//	cout <<"Bot1: " << Bot1 <<endl;
//	cout <<"Bot2: " << Bot2 <<endl;
//	cout <<"ADR2: " << ADR2 <<endl;
//	cout << "frac: " << (fac3*(Bot1*fac1*ADR1 + Bot2*fac2*ADR2 )) << endl;
	return fac3*(Bot1*fac1*ADR1 + Bot2*fac2*ADR2 );
}

std::complex<double> dblReggeMod::V12(double alp1, double alp2, double eta) const{

	if(alp1==alp2 ){return 0.0;}
	std::complex<double> res = CHGM(-alp1, 1.-alp1+alp2, -1/eta);
	res *= cgamma(alp1-alp2,0)/cgamma(-alp2,0);

//	cout << "alp1, alp2: " << alp1 <<", " <<alp2 <<endl;
//	cout << "res: " << res << endl;
	return res;
}

std::complex<double> dblReggeMod::DoubleRegge(int tau[2], double s, double si[2], double alp[2]) const{
	std::complex<double> ui (0,1);
	// signature factors:
	std::complex<double> x0  = 1/2.*((double)tau[0] + exp(-ui*M_PI*alp[0]));
	std::complex<double> x1  = 1/2.*((double)tau[1] + exp(-ui*M_PI*alp[1]));
	std::complex<double> x01 = 1/2.*((double)tau[0]*tau[1] + exp(-ui*M_PI*(alp[0]-alp[1])));
	std::complex<double> x10 = 1/2.*((double)tau[1]*tau[0] + exp(-ui*M_PI*(alp[1]-alp[0])));
	// double Regge vertices:
	double eta = S0*s/(si[0]*si[1]);
	std::complex<double> V0 = V12(alp[0], alp[1], eta);
	std::complex<double> V1 = V12(alp[1], alp[0], eta);
//	cout << "V0, V1: " << V0 <<", " <<V1 << endl;
//	cout << "x0, x1: " << x0 <<", "<<x1<< endl;
//	cout << "x01, x10: " <<x01<<", "<<x10<<endl;
		
	// combine pieces:
	std::complex<double> t1 = pow(s/S0,alp[1])*pow(si[0]/S0,alp[0]-alp[1])*x1*x01*V1;
	std::complex<double> t0 = pow(s/S0,alp[0])*pow(si[1]/S0,alp[1]-alp[0])*x0*x10*V0;

//	cout << "t0, t1: " << t0 <<", " <<t1 << endl;
	return (t0+t1)*cgamma(-alp[0],0)*cgamma(-alp[1],0);;
}


std::complex<double> dblReggeMod::cgamma(std::complex<double> z,int OPT) const{
	std::complex<double> ui (0,1);
	std::complex<double> g, infini= 1e308+ 0.0*ui; // z0,z1
	double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
	double na=0.0,t,x1 = 1,y1=0.0,sr,si;
	int j,k;


	static double a[] = {
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

	x = real(z);
	y = imag(z);

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

//	cout << "g: " << g << endl;
	return g;
}



double dblReggeMod::CHGM(double A, double B, double X) const{
	double A0=A, A1=A, X0=X, HG = 0.0;
	double TBA, TB, TA, Y0=0.0, Y1=0.0, RG, LA = (int) A, NL, R, M, INF = pow(10,300);
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
	if(HG != 0){return HG;}

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

//	cout << "HG: " << HG << endl;
	return HG;
}

void
dblReggeMod::updatePar( const AmpParameter& par ){

	// could do expensive calculations here on parameter updates
	//
}
