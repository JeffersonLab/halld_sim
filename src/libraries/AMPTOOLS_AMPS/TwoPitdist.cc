
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPitdist.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system) Elton 4/17/2017
// Version for sigma f0(500) production  Elton 9/9/2018
// Pealed off exponential dependence assuming it factorizes from the W dependence.

TwoPitdist::TwoPitdist( const vector< string >& args ) :
UserAmplitude< TwoPitdist >( args )
{
  
  assert( args.size() == 4 );
  Bslope = AmpParameter( args[0] );
  Bgen = AmpParameter( args[1] );
  m_daughters = pair< string, string >( args[2], args[3] );    // specify indices of pions in event
  
  // need to register any free parameters so the framework knows about them
  registerParameter( Bslope );
  registerParameter( Bgen );
  
  // make sure the input variables look reasonable
  assert( ( Bgen >= 1 ) && ( Bslope >= Bgen ) );     // Make sure generated value is lower than actual.         
}

complex< GDouble >
TwoPitdist::calcAmplitude( GDouble** pKin ) const
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
  
  /*GDouble Wpipi  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();*/

  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // Recoil is particle 3
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);      

    complex<GDouble> Arel(sqrt(exp(Bslope*t)/exp(Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return sqrt(exp^Bt) 
  

    // cout << " TwoPitdist" << " Bslope=" << Bslope << " Bgen=" << Bgen << " t=" << t <<  " Re(Arel)=" << real(Arel) << " imag(Arel)=" << imag(Arel) << endl; 
  return( Arel );
}

void
TwoPitdist::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


