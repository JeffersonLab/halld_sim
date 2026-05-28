#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <complex>

#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "breakupMomentumComplex.h"
#include "AMPTOOLS_AMPS/PiPiSWaveAMPK.h"


// S-wave dipion mass parametrization amplitude based on:
//   K.L. Au, D. Morgan and M.R. Pennington," Meson dynamics beyond the quark model: Study of final-state interactions",
//   PRD 35 1987
//   with the modification by I. Kachaev from [arxiv:2305.11711] 
//   Usage:
//   amplitude <reaction>::<sum>::<ampName> PiPiSWaveAMPK <id_daughter1> <id_daughter2> 

using namespace std;

PiPiSWaveAMPK::PiPiSWaveAMPK( const vector <string> &args ) : UserAmplitude<PiPiSWaveAMPK>(args)
{

   assert( args.size() == 2 );
   m_daughters = pair <string,string> (args[0],args[1]);

   setParametrization();

   // need to register a free parameters so the framework knows about it
   // registerParameter(m_mass);

}


void PiPiSWaveAMPK::setParametrization() 
{
  
  //Adler pole and residue in it
  s0 = -0.0074;  // AMP Table 1, M solution: s_0
  a11 = 0.1131; // AMP Table 1, M solution: f_2^2
  
  //polynom up to the 3rd order
  c11.resize(4);

  c11[0] = 0.0337;                 // AMP Table 1, M solution: c_11^0
  c11[1] = -0.3185;                // AMP Table 1, M solution: c_11^1
  c11[2] = -0.0942;                // AMP Table 1, M solution: c_11^2
  c11[3] = -0.5927;                // AMP Table 1, M solution: c_11^3

}



/////////////////////// Amplitude Calculation //////////////////////////

complex <GDouble> PiPiSWaveAMPK::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector pion1_P4, pion2_P4, dipion_P4, temp_P4;
  
  for( unsigned int ii = 0; ii < m_daughters.first.size(); ++ii ){

    char num[2]= {m_daughters.first[ii], '\0'};
    int index = atoi(num);

    temp_P4.SetPxPyPzE(pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0]);
    pion1_P4 += temp_P4;
    dipion_P4 += temp_P4;
  }
  
  for( unsigned int ii = 0; ii < m_daughters.second.size(); ++ii ){
    
    char num[2]= {m_daughters.second[ii], '\0'};
    int index = atoi(num);

    temp_P4.SetPxPyPzE(pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0]);
    pion2_P4 += temp_P4;
    dipion_P4 += temp_P4;
  }


  GDouble mass  = dipion_P4.M();
  GDouble s  = dipion_P4.M2();
   

  if (fabs(s - s0) < 1.e-6){
     mass += 1.e-6;
     s = mass*mass;
  }

  
  const complex <GDouble> qPiPi   = breakupMomentumComplex(mass, _piChargedMass,   _piChargedMass );
  const complex <GDouble> qPi0Pi0 = breakupMomentumComplex(mass, _piNeutralMass,   _piNeutralMass );
  const GDouble scale = (s / (4*_kaonMeanMass*_kaonMeanMass)) - 1;


  complex <GDouble> rho11 = (qPiPi + qPi0Pi0) / mass;

  complex <GDouble> imag(0.,1.);  
  complex <GDouble> M11(0.,0.),T11(0.,0.);

  //Add the Adler pole term
   M11 += a11/(s - s0);

  //Add the polynomial terms  
  for (unsigned int ii = 0; ii < c11.size(); ++ii) 
    M11 += c11[ii]*pow(scale, (int)ii);

  //Now calculate the T_{11} matrix element using the T11 to M11 relation
  T11 = 1./(M11 - imag*rho11); 


  return T11;
}




//This function may be used instead of the separate amplitude ' breakupMomentumComplex' 
template<typename mType> complex<GDouble> PiPiSWaveAMPK::breakupMom(mType m, GDouble mDec1, GDouble mDec2) const{
   complex<GDouble> result = PiPiSWaveAMPK::phaseSpaceFac(m, mDec1, mDec2)*m/GDouble(2.);
   return result;  
}


#ifdef GPU_ACCELERATION
void PiPiSWaveAMPK::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUPiPiSWaveAMPK_exec( dimGrid,  dimBlock, GPU_AMP_ARGS);
  
}
#endif 
