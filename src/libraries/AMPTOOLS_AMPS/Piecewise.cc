

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Piecewise.h"

Piecewise::Piecewise( const vector< string >& args ) :
UserAmplitude< Piecewise >( args )
{
  
  assert( args.size() == ( (2*atoi(args[2].c_str()) ) +6) );

  m_massMin   = atof( args[0].c_str() );
  m_massMax   = atof( args[1].c_str() );
  m_nBins     = atoi( args[2].c_str() );
  m_daughter1 = atoi( args[3].c_str() );
  m_daughter2 = atoi( args[4].c_str() );

  m_suffix = args[5];

  one = complex<GDouble>(1,0);
  zero = complex<GDouble>(0,0);

  m_width = (double)(m_massMax-m_massMin)/m_nBins;

  for(int i=0; i<m_nBins; i++) {
     string nameRe = "pcwsBin_" + to_string(i) + "Re" + m_suffix;
     string nameIm = "pcwsBin_" + to_string(i) + "Im" + m_suffix;

	// USING std::vectors
     m_paramsRe.push_back( AmpParameter( args[(2*i)+6] ) );
     m_paramsRe[i].setName( nameRe );
     m_paramsIm.push_back( AmpParameter( args[(2*i+1)+6] ) );
     m_paramsIm[i].setName( nameIm );
//     registerParameter( m_paramsRe[i] );
//     registerParameter( m_paramsIm[i] );

	// USING arrays
//     m_paramsRe[i].setName( nameRe );
//     m_paramsRe[i] = AmpParameter( args[(2*i)+5] );
//     m_paramsIm[i].setName( nameIm );
//     m_paramsIm[i] = AmpParameter( args[(2*i+1)+5] );
//     registerParameter( m_paramsRe[i] );
//     registerParameter( m_paramsIm[i] );
     
//     cout << m_params.back().name() << endl;
//     cout << m_params[i] << endl;
  }
  for(int i=0; i<m_nBins; i++) {
     registerParameter( m_paramsRe[i] );
     registerParameter( m_paramsIm[i] );
  }
}

complex< GDouble >
Piecewise::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2;
  
  P1.SetPxPyPzE( pKin[m_daughter1][1], pKin[m_daughter1][2],
                      pKin[m_daughter1][3], pKin[m_daughter1][0] );
  P2.SetPxPyPzE( pKin[m_daughter2][1], pKin[m_daughter2][2],
                      pKin[m_daughter2][3], pKin[m_daughter2][0] );

  GDouble mass = (P1+P2).M();

  int tempBin = 0;

  for(int i=0; i<m_nBins; i++) {
    if(mass>(m_massMin+(i*m_width)) && mass<(m_massMin+((i+1)*m_width)))
       tempBin = i;
  }
  
  return (complex<double>(m_paramsRe[tempBin],m_paramsIm[tempBin]));
 
}

void
Piecewise::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

void Piecewise::init(){
//  for(int i=0; i<m_nBins; i++) {
//     registerParameter( m_paramsRe[i] );
//     registerParameter( m_paramsIm[i] );
//  }
}

//#ifdef GPU_ACCELERATION
//void
//Piecewise::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
//  
//  // use integers to endcode the string of daughters -- one index in each
//  // decimal place
//  
//  int daught1 = atoi( m_daughters.first.c_str() );
//  int daught2 = atoi( m_daughters.second.c_str() );
//  
//  GPUBreitWigner_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, 
//                       m_mass0, m_width0, m_orbitL, daught1, daught2 );
//
//}
//#endif //GPU_ACCELERATION

