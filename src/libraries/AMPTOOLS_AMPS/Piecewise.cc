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
  
  assert( args.size() == ( (2*atoi(args[2].c_str()) ) + 7 ) );

  m_massMin   = atof( args[0].c_str() );
  m_massMax   = atof( args[1].c_str() );
  m_nBins     = atoi( args[2].c_str() );
  m_daughter1 = atoi( args[3].c_str() );
  m_daughter2 = atoi( args[4].c_str() );

  m_suffix = args[5]; // in case more than one piecewise amplitude is used in the cfg file, this string may contain a suffix to be added to all parameter names

  // switch between representation of complex parameters in Re/Im format and Mag/Phi format
  if(args[6] == "ReIm")
    m_represReIm = true;
  else if(args[6] == "MagPhi")
     m_represReIm = false;
  else
     cout << "ERROR: '" << args[6] << "' is not a defined mode for the Piecewise amplitude! Choose 'ReIm' or 'MagPhi'."
  
  m_width = (double)(m_massMax-m_massMin)/m_nBins;

  for(int i=0; i<m_nBins; i++) {
     string name1, name2;

     if(m_represReIm) {
        name1 = "pcwsBin_" + to_string(i) + "Re" + m_suffix;
        name2 = "pcwsBin_" + to_string(i) + "Im" + m_suffix;
     } else {
        name1 = "pcwsBin_" + to_string(i) + "Mag" + m_suffix;
        name2 = "pcwsBin_" + to_string(i) + "Phi" + m_suffix;
     }



     m_params1.push_back( AmpParameter( args[(2*i)+7] ) );
     m_params1[i].setName( name1 );
     m_params2.push_back( AmpParameter( args[(2*i+1)+7] ) );
     m_params2[i].setName( name2 );
  }

  for(int i=0; i<m_nBins; i++) {
     registerParameter( m_params1[i] );
     registerParameter( m_params2[i] );
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
  
  return (complex<double>(m_params1[tempBin],m_params2[tempBin]));
 
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

