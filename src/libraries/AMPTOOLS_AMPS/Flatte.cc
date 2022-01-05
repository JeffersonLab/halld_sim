#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Flatte.h"

// Flatte-type amplitude based on:
//   S.M. Flatte," Coupled-channel analysis of the pi eta and KKbar systems near KKbar threshold"
//   PLB 63 1976
//
//   <g1> and <g2> are the couplings of a resonance with mass <mass> decaying to <id_daughter1> and <id_daughter2>. 
//   Masses of the daughter particles for both relevant channels must be supplied (<chan1mass1>, ...) and finally the 
//   channel to be used for calculation is defined by <channel>, which can take values 1 or 2.
//   
//   Usage:
//   amplitude <reaction>::<sum>::<ampName> Flatte <mass> <g1> <g2> <id_daughter1> <id_daughter2> <chan1mass1> <chan1mass2> <chan2mass1> <chan2mass2> <channel>

using namespace std;

Flatte::Flatte( const vector< string >& args ) :
   UserAmplitude< Flatte >(args)
{

   assert( args.size() == 10 );

   m_mass=AmpParameter(args[0]); 
   m_g1=AmpParameter(args[1]); 
   m_g2=AmpParameter(args[2]); 
   m_daughter1 = atoi(args[3].c_str());
   m_daughter2 = atoi(args[4].c_str());
   m_mass11 = atof(args[5].c_str());
   m_mass12 = atof(args[6].c_str());
   m_mass21 = atof(args[7].c_str());
   m_mass22 = atof(args[8].c_str());
   m_chan = atoi(args[9].c_str());


   // need to register any free parameters so the framework knows about them
   registerParameter(m_mass);
   registerParameter(m_g1);
   registerParameter(m_g2);

}


complex< GDouble > Flatte::calcAmplitude( GDouble** pKin, GDouble* userData ) const {
   TLorentzVector P1, P2;

   P1.SetPxPyPzE( pKin[m_daughter1][1], pKin[m_daughter1][2], pKin[m_daughter1][3], pKin[m_daughter1][0] );
   P2.SetPxPyPzE( pKin[m_daughter2][1], pKin[m_daughter2][2], pKin[m_daughter2][3], pKin[m_daughter2][0] );

   double curMass = (P1+P2).M();
   complex<double> imag(0.,1.);

   complex<double> gamma11 = (double)m_g1 * Flatte::breakupMom( curMass, m_mass11, m_mass12 );
   complex<double> gamma22 = (double)m_g2 * Flatte::breakupMom( curMass, m_mass21, m_mass22 );

   complex<double> gammaLow;
   if( (m_mass11+m_mass12) < (m_mass21+m_mass22) ) gammaLow = gamma11;
   else gammaLow = gamma22;

   complex<double> gamma_j;
   if(m_chan==1) gamma_j = gamma11;
   else if(m_chan==2) gamma_j = gamma22;
   else cout << "ERROR: possible channel indices for Flatte amplitude are 1 or 2!" << endl;

   complex<double>  result = (double)m_mass * sqrt( gammaLow*gamma_j ) / ( (double)m_mass*(double)m_mass - curMass*curMass - imag * (double)m_mass * (gamma11+gamma22) );
   return result;
}


//void Flatte::calcUserData( GDouble** pKin, GDouble* userData ) const {
//}


complex<double> Flatte::phaseSpaceFac(double m, double mDec1, double mDec2) const{

   complex<double> result(0.,0.);

   if(fabs(m) < 1e-8) {
      std::cout << "Mass " << m << " is too close to 0. Cant calculate phasespace factor: \n set mass to 1.*e-10"  << endl;
      m=1.e-10;
   }

   double termPlus=(mDec1+mDec2)/m;
   double termMinus=(mDec1-mDec2)/m;
   double tmpVal=(1.-termPlus*termPlus) * (1.-termMinus*termMinus);
   if(tmpVal>=0.) result = complex<double>(std::sqrt(tmpVal), 0.);
   else result = complex<double>(0., sqrt(-tmpVal));   
   return result;
}

template<typename mType>
complex<double> Flatte::breakupMom(mType m, double mDec1, double mDec2) const{
   complex<double> result = Flatte::phaseSpaceFac(m, mDec1, mDec2)*m/2.;
   return result;  
}


#ifdef GPU_ACCELERATION
void
Flatte::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

   Flatte_exec( dimGrid,  dimBlock, GPU_AMP_ARGS,
         m_mass, m_g1, m_g2, m_daughter1, m_daughter2);

}
#endif //GPU_ACCELERATION
