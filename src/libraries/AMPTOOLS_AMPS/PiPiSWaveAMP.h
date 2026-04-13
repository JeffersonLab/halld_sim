#if !defined(PIPISWAVEAMP)
#define PIPISWAVEAMP

#include "particleType.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/Amplitude.h"
#include "GPUManager/GPUCustomTypes.h"


#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <iostream>


#ifdef GPU_ACCELERATION
void PiPiSWaveAMP_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble m_mass, GDouble m_g1, GDouble m_g2, int m_daughter1, int m_daughter2 );
#endif // GPU_ACCELERATION



using std::complex;
using namespace std;

class Kinematics;

class PiPiSWaveAMP : public UserAmplitude<PiPiSWaveAMP>{

   public:

      PiPiSWaveAMP() : UserAmplitude <PiPiSWaveAMP>() { }
      PiPiSWaveAMP( const vector<string> &args );
      ~PiPiSWaveAMP(){}

      string name() const { return "PiPiSWaveAMP"; }
      void setParametrizationAMPK(); 
      complex <GDouble> calcAmplitude( GDouble** pKin ) const;

#ifdef GPU_ACCELERATION

     void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
     bool isGPUEnabled() const { return false; }

#endif // GPU_ACCELERATION

   private:

      pair< string, string > m_daughters;    

      int _vesSheet;
      std::vector <GDouble> _sP; 
      std::vector <GDouble> _a11; 
      std::vector <GDouble> _a12; 
      std::vector <GDouble> _a21; 
      std::vector <GDouble> _a22; 
      std::vector <GDouble> _c11; 
      std::vector <GDouble> _c12; 
      std::vector <GDouble> _c21; 
      std::vector <GDouble> _c22; 

  
      const GDouble _piChargedMass = ParticleMass(PiPlus); //0.13957039;
      const GDouble _piNeutralMass = ParticleMass(Pi0); //0.13497680;
      const GDouble _kaonChargedMass = ParticleMass(KPlus); //0.49367700;
      const GDouble _kaonNeutralMass = 0.5*(ParticleMass(KLong) + ParticleMass(KShort)); //0.49761400;
      const GDouble _kaonMeanMass = 0.5*(_kaonChargedMass + _kaonNeutralMass); //0.5*(0.49367700 + 0.49761400);


  
      complex<GDouble> phaseSpaceFac( GDouble m, GDouble mDec1, GDouble mDec2 ) const;
      template<typename mType>
      complex<GDouble> breakupMom( mType m, GDouble mDec1, GDouble mDec2 ) const;

};

#endif
