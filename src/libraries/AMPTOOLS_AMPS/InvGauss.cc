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
#include "AMPTOOLS_AMPS/InvGauss.h"


// Class for modeling a t-distribution with inverse Gaussian
// three arguments are parsed:
// [0] <Lambda> double-valued, positive
// [1] <Mu> double-valued, positive
// [2] <m_daughters> are indices of the baryon decay products written in one string, e.g. "134"
 

InvGauss::InvGauss( const vector< string >& args ) : UserAmplitude< InvGauss >( args )
{  
  assert( args.size() == 3 );

  Lambda = AmpParameter( args[0] );
  Mu = AmpParameter( args[1] );
  m_daughters =  args[2];

  
  // register any free parameters so the framework knows about them
  registerParameter( Lambda );
  registerParameter( Mu );
  
  // make sure that all the input variables are positive
  assert( ( Lambda > 0) && ( Mu > 0 ) && m_daughters.length() > 0);     
}




complex<GDouble> InvGauss::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector recoil_baryon, tmp;
  TLorentzVector target( 0, 0, 0, ParticleMass(Proton) );
  GDouble minust, inverseGaussian;

    
  //identify four-momentum of the recoiled baryon and calculate a momentum transfer
  
  for( unsigned int ii = 0; ii < m_daughters.size(); ++ii ){

    char num[2]= {m_daughters[ii], '\0'};
    int index = atoi(num);

    tmp.SetPxPyPzE(pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0]);
    recoil_baryon += tmp;
  }
  
  
  minust = (recoil_baryon - target).M2()*(-1.);

  inverseGaussian = sqrt(Lambda/(2*TMath::Pi()*pow(minust,3)))*exp((-Lambda*(minust-Mu)*(minust-Mu))/(2.*Mu*Mu*minust));
  
  complex<GDouble> AmpRe(sqrt(inverseGaussian), 0.);  

  //  cout << " -t = " << minust <<  "\t Re(AmpRe)=" << real(AmpRe) << "\t Imag(AmpRe)=" << imag(AmpRe) << endl;
  
  return( AmpRe );
}

void InvGauss::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


