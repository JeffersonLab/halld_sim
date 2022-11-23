
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiEtas_tdist.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system) Elton 4/17/2017
// Version for sigma f0(500) production  Elton 9/9/2018
// Pealed off exponential dependence assuming it factorizes from the W dependence.

TwoPiEtas_tdist::TwoPiEtas_tdist( const vector< string >& args ) :
UserAmplitude< TwoPiEtas_tdist >( args )
{
  
  assert( args.size() == 5 );
  ThetaSigma = AmpParameter( args[0] );
  Phase = AmpParameter( args[1] );    // convert phase to radians
  Bgen = AmpParameter( args[2] );
  m_daughters = pair< string, string >( args[3], args[4] );    // specify indices of pions in event
  
  // need to register any free parameters so the framework knows about them
  registerParameter( ThetaSigma );
  registerParameter( Phase );
  registerParameter( Bgen );
  
  // make sure the input variables look reasonable
  assert( ( ThetaSigma >= 0) &&( Bgen >= 2 ) && (Phase >=0 && Phase <=180) );     // Make sure generated value is lower than actual.         
}

complex< GDouble >
TwoPiEtas_tdist::calcAmplitude( GDouble** pKin ) const
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


  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // Recoil is particle 3
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);  
  GDouble Eg = Ptot.E();
  GDouble tpar = (mass1*mass1/(2*Eg)) * (mass1*mass1/(2*Eg));

  GDouble Thpipi = -t > tpar? (180/PI)*G_SQRT( (-t-tpar)/(Eg*Ppipi) ): 0;

  // complex<GDouble> Arel(G_SQRT(exp(Bslope*t)/exp(Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return G_SQRT(exp^Bt) 
  //  complex<GDouble> Arel(G_SQRT(-t*exp(Bslope*t)/exp(Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return G_SQRT(-t*exp^Bt)   Add -t factor for pions 
  
  // Estimate of k2 sinthe / (-t) * F_strong(-t)  . PRC 80 055201 (2009) Eq. 4 and Fig 6.

  // complex<GDouble> Arel( (Eg*Eg/(-t)) * Thpipi * G_SQRT(exp(-Thpipi*Thpipi/(2*0.45*0.45))) /exp(Bgen*t),0. );   // 1/-t for photon exchange.

  complex<GDouble> Arel;
  complex<GDouble> ImagOne(0,1);


    Double_t arg = Thpipi*Thpipi/(2*ThetaSigma*ThetaSigma)<100? Thpipi*Thpipi/(2*ThetaSigma*ThetaSigma) : 0;
    Arel = arg > 0? G_SQRT(Thpipi) * G_SQRT(exp(-arg) /exp(Bgen*t)) :0;   // Can Change phase of sigma relative to Primakoff
    // cout << " arg=" << arg << " num=" << Thpipi * G_SQRT(exp(-arg)) << " den=" << exp(Bgen*t) << " Re(ImagOne)=" << real(ImagOne) << " imag(ImagOne)=" << imag(ImagOne)  << endl;
  // Adjust phase for this amplitude relative to other amplitudes
    Arel = Arel * (G_COS(Phase*PI/180.) + ImagOne*G_SIN(Phase*PI/180.));

    // cout << "TwoPiEtas_tdist" << " ThetaSigma=" << ThetaSigma << " Phase=" << Phase << " Bgen=" << Bgen << " t=" << t <<  " Re(Arel)=" << real(Arel) << " imag(Arel)=" << imag(Arel)  << " Eg=" << Eg << " tpar=" << tpar << " Wpipi=" << Wpipi << " Ppipi=" << Ppipi << " Thpipi=" << Thpipi << endl; 
  return( Arel );
}

void
TwoPiEtas_tdist::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


