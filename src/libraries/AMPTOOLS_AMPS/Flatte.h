#if !defined(FLATTE)
#define FLATTE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <iostream>

#ifdef GPU_ACCELERATION

void Flatte_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
      GDouble m_mass, GDouble m_g1, GDouble m_g2,
      int m_daughter1, int m_daughter2 );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class Flatte : public UserAmplitude< Flatte >{

   public:

      Flatte() : UserAmplitude< Flatte >() { }

      Flatte( const vector< string >& args );

      ~Flatte(){}

      string name() const { return "Flatte"; }

      complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userData ) const;

#ifdef GPU_ACCELERATION

      void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

      bool isGPUEnabled() const { return false; }

#endif // GPU_ACCELERATION

   private:

      AmpParameter m_mass;	
      AmpParameter m_g1;	
      AmpParameter m_g2;	
      int m_daughter1;
      int m_daughter2;  
      double m_mass11;
      double m_mass12;
      double m_mass21;
      double m_mass22;
      int m_chan;

      complex<double> phaseSpaceFac( double m, double mDec1, double mDec2 ) const;
      template<typename mType>
         complex<double> breakupMom( mType m, double mDec1, double mDec2 ) const;

};

#endif
