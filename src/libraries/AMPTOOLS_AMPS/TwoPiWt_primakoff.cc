
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"
#include "particleType.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiWt_primakoff.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system)
// Elton 4/17/2017

Double_t sigma_ggpi0pi0_func (Double_t *x, Double_t *par){
    
    // Parameterization for data from Masiske Phys Rev D 41 (1990) 3324.
    // Returns cross section in units of nb
    
    // constants
    // Double_t const PI = 3.14159;
    // parin[0] = A;              // parameter 1: normalization
    // parin[x] = mu;                // parameter 2: mean of Fermi Function. simplify to 2 parameters
    // parin[1] = kT              // parameter 3: transition parameter
    
  Double_t MPI = ParticleMass(Pi0);
    // Double_t PI=3.14159;
    
    Double_t A  = par[0];
    Double_t mu = 2*MPI;
    Double_t kT = par[1];
    Double_t Wpipi = x[0] ;
    Double_t f;
    
    if (Wpipi < 2*MPI) {
        f = 0.5;
    }
    else {
        f = 1 / (exp((Wpipi-mu)/kT) + 1);	                 // Use Fermi function. Data approx flat out to about 0.8 GeV, ignore beyond that. Data for costhe<0.8
    }
    
    
    // cout << " Wpipi=" << Wpipi << " A=" << A << " mu=" << mu << " kT=" << kT << " f=" << f << endl;
    
    return	2*A*(0.5-f)/0.8;                                     // Convert to nb, assuming isotropic angular distribution of pions. Correct for 80% coverage
}

Double_t sigma_ggpipi_func (Double_t *x, Double_t *par){
    
    // Parameterization from Rory for cross section for gamma gamma -> pi+ pi-
    // Returns cross section in units of nb/sr
    
    // constants
    // Double_t const PI = 3.14159;
  // Double_t MPI = ParticleMass(Pi0);     // pi0
  Double_t MPI = ParticleMass(PiPlus);  // pi+
    Double_t W0 = 0.3;
    
    Double_t expon = par[0];
    // Double_t par2 = par[1];
    Double_t Wpipi = x[0] ;
    Double_t f;
    
    if (Wpipi < 2*MPI) {
        f = 0;
    }
    else if (Wpipi < W0) {
        f = 300./(0.6*4.*PI)*(Wpipi-2.*MPI)/(W0-2.*MPI);           // linear rise, isotropic, CB data only 60% coverage
    }
    else {
        f = 300./(0.6*4.*PI)*pow(W0/Wpipi,expon);	                 // power fall off, isotropic
    }
    
    return	f;
}


Double_t ff_func (Double_t *x, Double_t *par){
    
    // return the nuclear form factor accourding to 2 parameter Fermi distribution
    // See Journall of Research of the National Bureau of Standards - B Mathenatics and Mathematical Physics
    // Vol. 70B, No. 1, Jan-Mar 1966. "The Form Factor of the Fermi Model Spatial Distribution," by Maximon and Schrack
    //
    // Function is a function of q, the three-momentum transfer to the nucleus.
    // Input argument is t
    // Note that q is the 3-vector momentum, but for low -t, q ~ G_SQRT(-t).
    
    // constants
    // Double_t alpha = 1/137.;
    // Double_t pi = 3.14159;
    Double_t hbarc = 0.19733;                  // GeV*fm
    Double_t q = G_SQRT(x[0])/hbarc;          // take q to be  in fm^-1. Input variable is positive (-t)
    
    Double_t R0  = par[0];  				 // half-density radius
    Double_t a0 = par[1];                    // skin or diffuseness parameter
    
    Double_t rho0;
    Double_t sum=0;
    Int_t jmax=4;
    for (Int_t j=1;j<jmax;j++) {                    // truncate after 3 terms, first term dominates.
        Double_t sum1 =pow(-1.,j-1)*exp(-j*R0/a0)/(j*j*j);
        sum += sum1;
        // cout << "jmax=" << jmax << " j=" << j << " R0=" << R0 << " a0=" << a0 << " sum1=" << sum1 << " sum=" << sum << endl;
    }
    
    rho0 = (4*PI*R0/3)*( PI*a0*PI*a0  + R0*R0 + 8*PI*a0*a0*a0*sum);
    rho0 = 1/rho0;
    
    Double_t f = 0;
    
    f = (4*PI*PI*rho0*a0*a0*a0)/(q*q*a0*a0*sinh(PI*q*a0)*sinh(PI*q*a0))
    	* (PI*q*a0 * cosh(PI*q*a0) * G_SIN(q*R0) - q*R0 *G_COS(q*R0) * sinh(PI*q*a0) )
    	+ 8*PI*rho0*a0*a0*a0*sum;
    
    // cout << " q=" << q << " f=" << f << endl;
    return f;
    
}


Double_t sigmat_func (Double_t *x, Double_t *par){
    
    // return the cross section for Primakoff production of pi+pi-. CPP proposal PR12-13-008 Eq 8.
    // independent variable is momentum transfer -t, in GeV2;
    
    // constants
    Double_t alpha = 1/137.;
    // Double_t pi = 3.14159;
    Double_t betapipi = 0.999;      // beta=0.999 (W=0.3), beta=0.997 (W=0.4), beta= 0.986 (W=1.0 GeV)
    // Double_t sigmagg = 1;         // take sigma (gg -> pi+ pi-) = 1 for now
    // Double_t Z = 50;              // Z of Sn, target
    Double_t Z = 82;              // Z of Pb, target
    Double_t coef = 4*alpha*Z*Z/(PI);
    
    Double_t Wpipi = par[0];
    Double_t Eg = par[1];
    Double_t t = -x[0] ;     // theta of pipi system in the lab. Input Variable is positive (i.e. -t)
    
    Double_t xin[2];
    Double_t parin[2];
    
    xin[0] = -t;                       // input variable to ff is positive (-t)
    parin[0] = par[2];  	       // half-density radius, fm
    parin[1] = par[3];                 // diffuseness paramter, fm
    
    Double_t FF = ff_func(xin,parin);
    
    // include other masses here
    
    Double_t m1 = 0;      // mass of photon, incident beam
    // Double_t m2 = 108.*0.931494;    // mass of 116Sn, target
    Double_t m2 = 208.*0.931494;    // use Pb mass because it is in the particle list
    Double_t m3 = Wpipi;   // mass of 2pi system, scattered system
    Double_t m4 = m2;     // recoil target
    
    
    Double_t f = 0;
    
    Double_t s = m2*m2 + 2*Eg*m2;
    Double_t sqrts = s > 0? G_SQRT(s) : 0;
    if (s < 0) {
        cout << "*** sigma_func: s =" << s << " < 0!" << endl;
        return f;
    }
   
    Double_t E1cm = (s+m1*m1-m2*m2)/(2*sqrts);
    // Double_t E2cm = (s+m2*m2-m1*m1)/(2*sqrts);
    Double_t E3cm = (s+m3*m3-m4*m4)/(2*sqrts);
    // Double_t E4cm = (s+m4*m4-m3*m3)/(2*sqrts);
    
    Double_t p1cm = E1cm*E1cm - m1*m1? G_SQRT(E1cm*E1cm - m1*m1) : 0;
    // Double_t p2cm = E2cm*E2cm - m2*m2? G_SQRT(E2cm*E2cm - m2*m2) : 0;
    Double_t p3cm = E3cm*E3cm - m3*m3? G_SQRT(E3cm*E3cm - m3*m3) : 0;
    // Double_t p4cm = E4cm*E4cm - m4*m4? G_SQRT(E4cm*E4cm - m4*m4) : 0;
    
    Double_t arg = (m1*m1-m3*m3-m2*m2+m4*m4)/(2*sqrts);
    Double_t t0 = arg*arg - (p1cm - p3cm)*(p1cm - p3cm);
    // Double_t t1 = arg*arg - (p1cm + p3cm)*(p1cm + p3cm);
    
    Double_t betastar = Eg/(Eg + m2);
    Double_t gammastar = (Eg + m2)/sqrts;
    Double_t betapipicm = p3cm/E3cm;
    
    // Double_t thepipicm = t0 -t > 0? G_SQRT((t0 -t)/(p1cm*p3cm)) : 0;
    
    Double_t conv = 1./(gammastar*(1 + betastar/betapipicm));
    
    if (-t > -t0) {
        f = (coef/2)* Eg*Eg*Eg*Eg * (t0-t)* betapipi*betapipi * (FF*FF/(t*t))*conv*conv*conv*conv/(p1cm*p3cm*p1cm*p3cm);
    }
    else   {
        f = 0;
    }
    
    // if (f <= 0) cout << " t=" << t << " t0=" << t0 << " betastar=" << betastar << " gammastar=" << gammastar << " betapipicm=" 
    //   << betapipicm << " FF=" << FF << " f=" << f << endl;
    return	f;
}


TwoPiWt_primakoff::TwoPiWt_primakoff( const vector< string >& args ) :
UserAmplitude< TwoPiWt_primakoff >( args )
{
  
  assert( args.size() == 6);
	m_par1 = AmpParameter( args[0] );
	m_par2 = AmpParameter( args[1] );
	Bgen = AmpParameter( args[2] );
        mtmax = AmpParameter( args[3] );
	m_daughters = pair< string, string >( args[4], args[5] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_par1 );
  registerParameter( m_par2 );
  registerParameter( Bgen );
  registerParameter( mtmax );
  
  // make sure the input variables look reasonable
  // assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );  
  assert( Bgen >= 1);   
  assert( mtmax > 0);         
}

complex< GDouble >
TwoPiWt_primakoff::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp, Precoil;
  
  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
    
    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );    // pi+ is index 1
    P1 += Ptemp;
    Ptot += Ptemp;

    /* cout << " 1i=" << i << " num=" << num << " index=" << index  << " P1.M=" << P1.M() << endl;
    P1.Print();
    Ptot.Print();*/
  }
  
  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
    
    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );    // pi- is index 2
    P2 += Ptemp;
    Ptot += Ptemp;   

    /* cout << " 2i=" << i << " num=" << num << " index=" << index << " P2.M=" << P2.M() << endl;
    P2.Print();
    Ptot.Print();*/
  }
  
  GDouble Wpipi  = Ptot.M();
  GDouble mass1 = P1.M();
  // GDouble mass2 = P2.M();
  GDouble Ppipi = Ptot.E() > Wpipi? G_SQRT(Ptot.E()*Ptot.E() - Wpipi*Wpipi): 0;
  GDouble Thetapipi = Ptot.Theta()*180./PI;

  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // Recoil is particle 3
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);      

  // cout << "Precoil.M()=" << Precoil.M() << " T=" << Precoil.E() - Precoil.M() << " t=" << t << endl; Precoil.Print();cout << endl << endl;

  // call sigma (gamma gamma -> pi pi) cross section

  
    Int_t const npar = 4;
    Double_t xin[1];
    xin[0] = Wpipi;                // W, 2pi mass
    Double_t Eg = pKin[0][0];          // incident photon energy
    Double_t parin[npar];
    // parin[0] A= 9.6;      // pi0 fit to data for costhe<0.8  Marsiske Phys Rev D 41 (1990) 3324
    // parin[x]  mu= 2*MPI;     // pi0 should go to zero at threshold. Move to function
    // parin[1] kT= 0.028;      // pi0 transition over about 100 MeV.

    // parin[0] = 1.29;              // charged parameter 1: exponent
    // parin[1] = 0.;                // charged parameter 2: par2 (spare)
    parin[0] = m_par1;              // parameter 1: exponent
    parin[1] = m_par2;                // parameter 2: par2 (spare)
    // Double_t Wmin=0.2 ;
    // Double_t Wmax=0.8;


    // cout << " TwoPiWt_primakoff: m_par1=" << m_par1 << " m_par2=" << m_par2 << " Bgen=" << Bgen << endl;

    GDouble sig_ggpipi = sigma_ggpi0pi0_func(xin,parin);

    parin[0] = Wpipi;
    parin[1] = Eg;
    // Double_t R0 = 6.62;   // Pb half-density radius, fm
    // Double_t a0 = 0.546;   // Pb difuseness parameter, fm
    // Double_t R0  = 5.358;   // Sn half-density radius, fm
    // Double_t a0 = 0.550;   // Sn difuseness parameter, fm
    parin[2] = 6.62;   
    parin[3] = 0.546;  
    xin[0] = -t;             // input positive value of t
    
    GDouble xnorm = 0.001;
    GDouble sigmat = sigmat_func (xin,parin) * xnorm;    // normlize amplitude to about unity

  GDouble tpar = (mass1*mass1/(2*Eg)) * (mass1*mass1/(2*Eg));
  GDouble Thpipi = -t > tpar? (180/PI)*G_SQRT( (-t-tpar)/(Eg*Ppipi) ): 0;
  
  GDouble epsilon = 1e-7;
  complex<GDouble> RealOne(1,0);
  complex<GDouble> ImagOne(0,1);
  complex<GDouble> Csig;

  double arg = Bgen*t > -100? Bgen*t : -100;   // limit exponential
  Csig = isfinite(sigmat*sig_ggpipi/Wpipi/exp(arg))? G_SQRT(sigmat*sig_ggpipi/Wpipi/exp(arg)) * RealOne : 0;    // Return complex double, sqrt (cross section). Divide out exponential
  if (-t > mtmax) Csig = 0;      // eliminate events at high t with large weights


  // cout << "calcAmplitude: 2pi mass=" << Wpipi << " Eg=" << Eg << " t=" << t << " BGen=" << Bgen << " exp(Bgen*t)=" << exp(Bgen*t) << " m_par1=" << m_par1 << " m_par2=" 
  //<< m_par2 << " sig_ggpipi=" << sig_ggpipi << " sigmat=" << sigmat << " Thpipi=" << Thpipi << " Thetapipi=" << Thetapipi << " Csig=" << Csig << endl;

    return( Csig + epsilon);   // return non-zero value to protect from likelihood calculation
}

void
TwoPiWt_primakoff::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


