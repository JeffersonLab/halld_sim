#if !defined(PIPISWAVEAMPK)
#define PIPISWAVEAMPK

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


using std::complex;
using namespace std;


#ifdef GPU_ACCELERATION
void   GPUPiPiSWaveAMPK_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO );
#endif 



class Kinematics;

class PiPiSWaveAMPK : public UserAmplitude<PiPiSWaveAMPK>{

   public:

      PiPiSWaveAMPK() : UserAmplitude <PiPiSWaveAMPK>() { }
      PiPiSWaveAMPK( const vector<string> &args );
      ~PiPiSWaveAMPK(){}
      
      enum UserVars { uv_ampRe = 0, uv_ampIm = 0, kNumUserVars };
      unsigned int numUserVars() const { return kNumUserVars; }
      void calcUserVars( GDouble **pKin, GDouble *userVars ) const;  
  
      string name() const { return "PiPiSWaveAMPK"; }   
      void setParametrization(); 
      complex <GDouble> calcAmplitude( GDouble **pKin, GDouble *userVars ) const;

  
#ifdef GPU_ACCELERATION

     void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
     bool isGPUEnabled() const { return true; }

#endif 

   private:

      pair< string, string > m_daughters;    

      GDouble s0; 
      GDouble a11; 
      std::vector <GDouble> c11; 

  
      const GDouble _piChargedMass = ParticleMass(PiPlus); //0.13957039;
      const GDouble _piNeutralMass = ParticleMass(Pi0); //0.13497680;
      const GDouble _kaonChargedMass = ParticleMass(KPlus); //0.49367700;
      const GDouble _kaonNeutralMass = 0.5*(ParticleMass(KLong) + ParticleMass(KShort)); //0.49761400;
      const GDouble _kaonMeanMass = 0.5*(_kaonChargedMass + _kaonNeutralMass); //0.5*(0.49367700 + 0.49761400);

  
      template<typename mType>
      complex<GDouble> breakupMom( mType m, GDouble mDec1, GDouble mDec2 ) const;

};

#endif
