#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <complex>

#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "breakupMomentumComplex.h"
#include "AMPTOOLS_AMPS/PiPiSWaveAMP.h"


// S-wave dipion mass parametrization amplitude based on:
//   K.L. Au, D. Morgan and M.R. Pennington," Meson dynamics beyond the quark model: Study of final-state interactions",
//   PRD 35 1987
//   with the modification by I. Kachaev from [arxiv:2305.11711] 
//   Usage:
//   amplitude <reaction>::<sum>::<ampName> PiPiSWaveAMP <id_daughter1> <id_daughter2> 

using namespace std;

PiPiSWaveAMP::PiPiSWaveAMP( const vector <string> &args ) : UserAmplitude<PiPiSWaveAMP>(args)
{

   assert( args.size() == 2 );
   m_daughters = pair <string,string> (args[0],args[1]);
   
   // need to register a free parameters so the framework knows about it
   //   registerParameter(m_mass);

}


void PiPiSWaveAMP::setParametrizationAMPK() 
{
  _sP.resize(2);
  _a11.resize(2);
  _a12.resize(2);
  _a21.resize(2);
  _a22.resize(2);

  _c11.resize(5);
  _c12.resize(5);
  _c21.resize(5);
  _c22.resize(5);

  
  //pole #0  
  _sP[0] = -0.0074;  // AMP Table 1, M solution: s_0
  _a11[0] = 0.1131; // AMP Table 1, M solution: f_2^2
  _a22[0] = -0.3216; // AMP Table 1, M solution: f_2^3

  
  //pole #1  (eliminated)
  _sP[1] = 0.9828; // AMP Table 1, M solution: s_1
  _a11[1] = 0.;    // -0.1968*0.1968 for the sign-corrected value from AMP Table 1, M solution: f_1^1
  _a22[1] = 0.;    // -0.0154*0.0154 for the sign-corrected value from AMP Table 1, M solution: f_2^1


  //polynom: #0th order
  _c11[0] = 0.0337;                // AMP Table 1, M solution: c_11^0
  _c22[0] = 0.3010;                // AMP Table 1, M solution: c_22^0

  //polynom: #1st order
  _c11[1] = -0.3185;                // AMP Table 1, M solution: c_11^1
  _c22[1] = -0.5140;                // AMP Table 1, M solution: c_22^1

  //polynom: #2nd order
  _c11[2] = -0.0942;                // AMP Table 1, M solution: c_11^2
  _c22[2] = 0.1176;                // AMP Table 1, M solution: c_22^2

  //polynom: #3rd order
  _c11[3] = -0.5927;                // AMP Table 1, M solution: c_11^3
  _c22[3] = 0.5204;                // AMP Table 1, M solution: c_22^3

  //polynom: #4th order (eliminated)
  _c11[4] = 0.;                //  0.1957 from AMP Table 1, M solution: c_11^4
  _c22[4] = 0.;                // -0.3977 from AMP Table 1, M solution: c_22^4

  
  //All off-diagonal elements go to zero
  _a12 = {0.,0.,0.,0.,0.};
  _a21 = {0.,0.,0.,0.,0.};
  _c12 = {0.,0.,0.,0.,0.};
  _c21 = {0.,0.,0.,0.,0.};

}



/////////////////////// Amplitude Calculation //////////////////////////

complex <GDouble> PiPiSWaveAMP::calcAmplitude( GDouble** pKin ) const
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
   

  if (fabs(s - _sP.back()) < 1.e-6){
     mass += 1.e-6;
     s = mass*mass;
  }

  
  const complex <GDouble> qPiPi   = breakupMomentumComplex(mass, _piChargedMass,   _piChargedMass );
  const complex <GDouble> qPi0Pi0 = breakupMomentumComplex(mass, _piNeutralMass,   _piNeutralMass );
  const GDouble scale = (s / (4*_kaonMeanMass*_kaonMeanMass)) - 1;


  complex <GDouble> rho11 = (qPiPi + qPi0Pi0) / mass;

  complex <GDouble> imag(0.,1.);  
  complex <GDouble> M11(0.,0.),T11(0.,0.);

  //Sum over the poles
  for (unsigned int ii = 0; ii < _sP.size(); ++ii)
    M11 += _a11[ii]/(s - _sP[ii]);

  //Sum over the polynomial terms  
  for (unsigned int ii = 0; ii < _c11.size(); ++ii) 
    M11 += _c11[ii]*pow(scale, (int)ii);

  //Now calculate the T_{11} matrix element using only diagonal terms
  T11 = M11 - imag*rho11;
 
  complex <GDouble>  amp = 1./T11;

   return amp;
}




//This function may be used instead of the separate amplitude ' breakupMomentumComplex' 
template<typename mType> complex<GDouble> PiPiSWaveAMP::breakupMom(mType m, GDouble mDec1, GDouble mDec2) const{
   complex<GDouble> result = PiPiSWaveAMP::phaseSpaceFac(m, mDec1, mDec2)*m/GDouble(2.);
   return result;  
}


#ifdef GPU_ACCELERATION
void PiPiSWaveAMP::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  PiPiSWaveAMP_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, m_mass, m_g1, m_g2, m_daughter1, m_daughter2);

}
#endif //GPU_ACCELERATION
