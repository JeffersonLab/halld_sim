
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
#include "AMPTOOLS_AMPS/Lambda1520tdist.h"

// Class modeled after TwoPitdist amplitude function
// 
// Pealed off exponential dependence assuming it factorizes from the W dependence.

Lambda1520tdist::Lambda1520tdist( const vector< string >& args ) :
UserAmplitude< Lambda1520tdist >( args )
{
  
  assert( args.size() == 3 );
  Bslope = AmpParameter( args[0] );
  exponent = AmpParameter( args[1] );
  Bgen = AmpParameter( args[2] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( Bslope );
  registerParameter( Bgen );
  
  // make sure the input variables look reasonable
  assert( ( Bgen >= 1 ) && ( Bslope >= Bgen ) );     // Make sure generated value is lower than actual.         
}

complex< GDouble >
Lambda1520tdist::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector d1, d2;
  TLorentzVector target( 0, 0, 0, ParticleMass(Proton) );

  // get momentum transfer
  d1.SetPxPyPzE (pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);   // daughter1 is particle 2
  d2.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // daughter2 is particle 3
  GDouble t = (d1+d2-target).M2()*(-1.);

    complex<GDouble> Arel(sqrt(TMath::Power(t,exponent)*exp(-Bslope*t)/exp(-Bgen*t)),0.);  // Divide out generated exponential. This must be the same as in GammaZToXYZ.cc. Return sqrt(exp^Bt) 
  

    // cout << " Lambda1520tdist" << " Bslope=" << Bslope << " Bgen=" << Bgen << " t=" << t <<  " Re(Arel)=" << real(Arel) << " imag(Arel)=" << imag(Arel) << endl; 
  return( Arel );
}

void
Lambda1520tdist::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


