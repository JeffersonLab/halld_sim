#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/EtaPb_tdist.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system) Elton 4/17/2017
// Version for sigma f0(500) production  Elton 9/9/2018
// Pealed off exponential dependence assuming it factorizes from the W dependence.

EtaPb_tdist::EtaPb_tdist( const vector< string >& args ) :
UserAmplitude< EtaPb_tdist >( args )
{
  
  assert( args.size() == 3 );            // assume eta is first and target is last
  ThetaSigma = AmpParameter( args[0] );
  Phase = AmpParameter( args[1] );    // convert phase to radians
  Bgen = AmpParameter( args[2] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( ThetaSigma );
  registerParameter( Phase );
  registerParameter( Bgen );
  
  // make sure the input variables look reasonable
  assert( ( ThetaSigma >= 0) &&( Bgen >= 2 ) && (Phase >=0 && Phase <=180) );     // Make sure generated value is lower than actual.         
}

complex< GDouble >
EtaPb_tdist::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector PEta, Precoil, Ptot, PGamma;

  PGamma.SetPxPyPzE (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]);   // Gamma is particle 0
  PEta.SetPxPyPzE (pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);   // Eta is particle 1
  Precoil.SetPxPyPzE (pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);   // Recoil is particle 2
  Ptot = PEta + Precoil;

  /*cout  << "M=" << PEta.M() << " PEta="; PEta.Print();
  cout  << "M=" << Precoil.M() << " Precoil="; Precoil.Print();
  cout  << "M=" << Ptot.M() << " Ptot="; Ptot.Print(); cout << endl;*/


  // get momentum transfer
  GDouble Et = Precoil.E();
  GDouble Mt = Precoil.M();
  GDouble t = -2*Precoil.M()*(Et - Mt);  
  GDouble Eg = PGamma.E();
  GDouble MEta = PEta.M();
  GDouble tpar = (MEta*MEta/(2*Eg)) * (MEta*MEta/(2*Eg));
  GDouble peta = sqrt(PEta.E()*PEta.E() - MEta*MEta);

  GDouble ThEta = -t > tpar? (180/PI)*sqrt( (-t-tpar)/(Eg*peta) ): 0;   // assumes lab is also cm frame. 3% difference for Eg and peta in cm

  // complex<GDouble> Arel(sqrt(exp(Bslope*t)/exp(Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return sqrt(exp^Bt) 
  //  complex<GDouble> Arel(sqrt(-t*exp(Bslope*t)/exp(Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return sqrt(-t*exp^Bt)   Add -t factor for pions 
  
  // Estimate of k2 sinthe / (-t) * F_strong(-t)  . PRC 80 055201 (2009) Eq. 4 and Fig 6.

  // complex<GDouble> Arel( (Eg*Eg/(-t)) * ThEta * sqrt(exp(-ThEta*ThEta/(2*0.45*0.45))) /exp(Bgen*t),0. );   // 1/-t for photon exchange.

  complex<GDouble> Arel;
  complex<GDouble> ImagOne(0,1);


    Double_t arg = ThEta*ThEta/(2*ThetaSigma*ThetaSigma)<100? ThEta*ThEta/(2*ThetaSigma*ThetaSigma) : 0;
    Arel = arg > 0? ThEta * sqrt(exp(-arg) /exp(Bgen*t)) :0;   // Can Change phase of sigma relative to Primakoff
  // Adjust phase for this amplitude relative to other amplitudes
    Arel = Arel * (G_COS(Phase*PI/180.) + ImagOne*G_SIN(Phase*PI/180.));
    // cout << " t=" << t << " ThEta=" << ThEta << " arg=" << arg << " num=" << ThEta * sqrt(exp(-arg)) << " den=" << exp(Bgen*t) << " Re(Arel)=" << real(Arel) << " imag(Arel)=" << imag(Arel)  << endl;

    //cout << "EtaPb_tdist" << " ThetaSigma=" << ThetaSigma << " Phase=" << Phase << " Bgen=" << Bgen << " t=" << t 
    //     << " Re(Arel)=" << real(Arel) << " imag(Arel)=" << imag(Arel)  << " Eg=" << Eg << " tpar=" << tpar << " peta=" << peta << " ThEta=" << ThEta << endl; 

  return( Arel );
}

void
EtaPb_tdist::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


