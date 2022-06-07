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
  
  assert( args.size() == uint( (2*atoi(args[2].c_str()) ) + 6 ) );

  m_massMin   = atof( args[0].c_str() );
  m_massMax   = atof( args[1].c_str() );
  m_nBins     = atoi( args[2].c_str() );
  m_daughters = string( args[3] );

  m_suffix = args[4]; // in case more than one piecewise amplitude is used in the cfg file, this string may contain a suffix to be added to all parameter names

  // switch between representation of complex parameters in Re/Im format and Mag/Phi format
  if(args[5] == "ReIm")
    m_represReIm = true;
  else if(args[5] == "MagPhi")
     m_represReIm = false;
  else
     cout << "ERROR: '" << args[5] << "' is not a defined mode for the Piecewise amplitude! Choose 'ReIm' or 'MagPhi'." << endl;
  
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

     m_params1.push_back( AmpParameter( args[(2*i)+6] ) );
     m_params1[i].setName( name1 );
     m_params2.push_back( AmpParameter( args[(2*i+1)+6] ) );
     m_params2[i].setName( name2 );
  }

  for(int i=0; i<m_nBins; i++) {
     registerParameter( m_params1[i] );
     registerParameter( m_params2[i] );
  }
}

complex< GDouble >
Piecewise::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
	// convert double back to long for array index
#ifdef AMPTOOLS_GDOUBLE_FP64
  long* tempBin = (long*)&(userVars[uv_imassbin]);
#else
  int* tempBin = (int*)&(userVars[uv_imassbin]);
#endif
  
	complex<GDouble> ans(m_params1[*tempBin],m_params2[*tempBin]);
	if(!m_represReIm)
		ans = polar(fabs(GDouble(m_params1[*tempBin])),GDouble(m_params2[*tempBin]));
		
	return ans;
}

void
Piecewise::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector Ptemp, Ptot;
  
  for( unsigned int i = 0; i < m_daughters.size(); ++i ){
    string num; num += m_daughters[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    Ptot += Ptemp;
  }

  GDouble mass = Ptot.M();

  
#ifdef AMPTOOLS_GDOUBLE_FP64
  long tempBin = 0;
#else
  int tempBin = 0;
#endif
  
  for(int i=0; i<m_nBins; i++) {
    if(mass>(m_massMin+(i*m_width)) && mass<(m_massMin+((i+1)*m_width)))
       tempBin = i;
  }
  
  // from Matt: use the memory allocated to a double type user variable to write the bin index as a long int
  userVars[uv_imassbin] = *((GDouble*)&tempBin); 

}

void
Piecewise::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

#ifdef GPU_ACCELERATION
void
Piecewise::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

        // convert vector to array for GPU
        GDouble params1[m_nBins];
        GDouble params2[m_nBins];
        for(int i=0; i<m_nBins; i++){
                params1[i] = m_params1[i];
                params2[i] = m_params2[i];
        }
        GPUPiecewise_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, params1, params2, m_nBins, m_represReIm);

}
#endif //GPU_ACCELERATION
