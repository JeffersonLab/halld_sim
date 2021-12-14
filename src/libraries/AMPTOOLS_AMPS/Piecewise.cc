

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
  
  assert( args.size() == uint( (2*atoi(args[2].c_str()) ) +5) );

  m_massMin   = atof( args[0].c_str() );
  m_massMax   = atof( args[1].c_str() );
  m_nBins     = atoi( args[2].c_str() );
  m_daughters = string( args[3] );

  m_suffix = args[4];

  one = complex<GDouble>(1,0);
  zero = complex<GDouble>(0,0);

  m_width = (double)(m_massMax-m_massMin)/m_nBins;

  for(int i=0; i<m_nBins; i++) {
     string nameRe = "pcwsBin_" + to_string(i) + "Re" + m_suffix;
     string nameIm = "pcwsBin_" + to_string(i) + "Im" + m_suffix;

	// USING std::vectors
     m_paramsRe.push_back( AmpParameter( args[(2*i)+5] ) );
     m_paramsRe[i].setName( nameRe );
     m_paramsIm.push_back( AmpParameter( args[(2*i+1)+5] ) );
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

void
Piecewise::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector Ptot, Ptemp;
  
  for( unsigned int i = 0; i < m_daughters.size(); ++i ){
    
    string num; num += m_daughters[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    Ptot += Ptemp;
  }
  
  GDouble mass  = Ptot.M();

  int tempBin = 0;

  for(int i=0; i<m_nBins; i++) {
    if(mass>(m_massMin+(i*m_width)) && mass<(m_massMin+((i+1)*m_width)))
       tempBin = i;
  }
  
  userVars[uv_imassbin] = tempBin;
}

complex< GDouble >
Piecewise::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
	int tempBin = userVars[uv_imassbin];
	return (complex<double>(m_paramsRe[tempBin],m_paramsIm[tempBin]));
}

void
Piecewise::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

void Piecewise::init(){

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

