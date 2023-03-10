
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
#include "AMPTOOLS_AMPS/TwoPiW_brokenetas.h"

// Class modeled after BreitWigner amplitude function provided for examples with AmpTools.
// Dependence of swave 2pi cross section on W (mass of 2pi system) Elton 4/17/2017
// Version for sigma f0(500) production  Elton 9/9/2018
// Simple version for broken etas (3pi0 reconstructed as 2pi0) Elton 5/7/2020

TwoPiW_brokenetas::TwoPiW_brokenetas( const vector< string >& args ) :
UserAmplitude< TwoPiW_brokenetas >( args )
{
  
  assert( args.size() == 4 );
	m_par1 = AmpParameter( args[0] );
	m_par2 = AmpParameter( args[1] );
	m_daughters = pair< string, string >( args[2], args[3] );
  
  // need to register any free parameters so the framework knows about them
  // for brokenetas, parameters are Gmean and Gsigma of the Gaussian in GeV
  registerParameter( m_par1 );
  registerParameter( m_par2 );
  
  // make sure the input variables look reasonable
  // assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );
}

complex< GDouble >
TwoPiW_brokenetas::calcAmplitude( GDouble** pKin ) const
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
  
  //  Double_t Eg = pKin[0][0];          // incident photon energy
  GDouble Wpipi  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();

  // get momentum transfer
  Precoil.SetPxPyPzE (pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);   // Recoil is particle 3
  // next three lines commented out, unused variables
  //  GDouble Et = Precoil.E();
  //  GDouble Mt = Precoil.M();
  //  GDouble t = -2*Precoil.M()*(Et - Mt);  

  complex<GDouble> ImagOne(0,1);
    complex<GDouble> Aw;
    GDouble Gmean= m_par1;
    GDouble Gsigma=m_par2;
    GDouble ImPart=0;

    Aw = exp( -(Wpipi-Gmean)*(Wpipi-Gmean)/(2*Gsigma*Gsigma)) + ImPart*ImagOne;

    Aw = isfinite(real(Aw))  && isfinite(imag(Aw))? Aw : 0;   // protect against infitinites

    if (Wpipi < mass1+mass2) Aw = 0;
    
    // cout << "TwoPiW_brokenetas: calcAmplitude: 2pi mass=" << Wpipi << " Eg=" << Eg << " t=" << t << " Gmean=" << Gmean << " Gsigma=" << Gsigma << " AwNorm=" << std::norm(Aw) << " AwPhase=" << std::arg(Aw) << endl;
  
  return( Aw );
}

void
TwoPiW_brokenetas::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}


