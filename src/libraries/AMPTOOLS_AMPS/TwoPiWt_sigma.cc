
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "particleType.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiWt_sigma.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system) Elton 4/17/2017
// Version for sigma f0(500) production  Elton 9/9/2018

complex<double> Aw_func (Double_t *x, Double_t *par) {
  // Returns the amplitude for the sigma photon production. GlueX-doc-XXXX DRAFT October 1, 2018
  // "Amplitude generation of the $\sigma$ meson f$_0$(500)"
    
    double alpha1Re=par[0];
    double alpha1Im=par[1];
    double alpha2Re=par[2];
    double alpha2Im=par[3];
    double alpha3Re=par[4];
    double alpha3Im=par[5];
    double alpha4Re=par[6];
    double alpha4Im=par[7];
    double W=x[0];
    
    complex<double> Amp1, Amp2;
    complex<double> alpha1(alpha1Re,alpha1Im);
    complex<double> alpha2(alpha2Re,alpha2Im);
    complex<double> alpha3(alpha3Re,alpha3Im);
    complex<double> alpha4(alpha4Re,alpha4Im);
    
    
    // S-wave I=0 phase shift as a function of center of mass energy
    
    Double_t Mpi = ParticleMass(Pi0);
    Double_t Mk = ParticleMass(KPlus);
    Double_t s = W*W;            // center of mass energy ^2
    Double_t s0 = (2*Mk)*(2*Mk);
    
    if (s <= 4*Mpi*Mpi || s >= s0) return 0;
    
    Double_t k = s/4 - Mpi*Mpi >= 0? sqrt(s/4 - Mpi*Mpi) : 0;   // center-of-mass momentum
 
    // Reference Ananthanaryan, Phys Reports 353 (2001) 207, Appendix D.   tan(delta)
    // S-wave I=0 phase shift as a function of center of mass energy
    Double_t A00=0.225;
    Double_t B00=12.651;
    Double_t C00=-43.8454;
    Double_t D00=-87.1632;
    Double_t s00=0.715311;
    
    Double_t k2=k*k;
    Double_t f = (2*k/sqrt(s))* (A00 + B00*k2 + C00*k2*k2 + D00*k2*k2*k2)*(4*Mpi*Mpi-s00)/(s-s00);
    
    Double_t delta0 = atan(f);
    if (delta0 < 0) delta0 += 3.14159;

    // cout << " s=" << s << " delta0=" << delta0 << endl;
    
    Double_t sind = sin(delta0);
    Double_t cosd = cos(delta0);
    complex<double> expdelta(cos(delta0),sin(delta0)); 
 
    Amp1 = (W/(2*k)) * sind * expdelta * (alpha1 + alpha2*s);
    Amp2 = cosd * expdelta* (alpha3 + alpha4*s);
    complex<double> Amp = Amp1 + Amp2;
    
    return Amp;

}


TwoPiWt_sigma::TwoPiWt_sigma( const vector< string >& args ) :
UserAmplitude< TwoPiWt_sigma >( args )
{
  
  assert( args.size() == 4 );
	m_par1 = AmpParameter( args[0] );
	m_par2 = AmpParameter( args[1] );
	m_daughters = pair< string, string >( args[2], args[3] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_par1 );
  registerParameter( m_par2 );
  
  // make sure the input variables look reasonable
  // assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );
}

complex< GDouble >
TwoPiWt_sigma::calcAmplitude( GDouble** pKin ) const
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
  GDouble mass2 = P2.M();

  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // Recoil is particle 3
  // next three lines commented out, unused variables
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);  

  
    Int_t const npar = 8;
    Double_t xin[1];
    xin[0] = Wpipi;                // W, 2pi mass
    Double_t Eg = pKin[0][0];          // incident photon energy

    Double_t alpha1Re=0.378129;
    Double_t alpha1Im=0;
    Double_t alpha2Re=0.751557;
    Double_t alpha2Im=0;
    Double_t alpha3Re=0.244899;
    Double_t alpha3Im=0;
    Double_t alpha4Re=0.179788;
    Double_t alpha4Im=0;

    Double_t parin[npar];
    parin[0] = alpha1Re;
    parin[1] = alpha1Im;
    parin[2] = alpha2Re;
    parin[3] = alpha2Im;
    parin[4] = alpha3Re;
    parin[5] = alpha3Im;
    parin[6] = alpha4Re;
    parin[7] = alpha4Im;

    complex<GDouble> Aw = Aw_func(xin,parin);

    /*Double_t Bslope=3.7; 
    Double_t Bgen=6.0;

    Aw = Aw * exp(Bslope*t)/exp(6.0*t);         // Divide out generated exponential. Move to separate amplitude.*/ 

    if (Wpipi < mass1+mass2) Aw = 0;
    
    // cout << "TwoPiWt_sigma: calcAmplitude: 2pi mass=" << Wpipi << " Eg=" << Eg << " t=" << t << " AwNorm=" << std::norm(Aw) << " AwPhase=" << std::arg(Aw) << endl;
  
  return( Aw );
}

void
TwoPiWt_sigma::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


