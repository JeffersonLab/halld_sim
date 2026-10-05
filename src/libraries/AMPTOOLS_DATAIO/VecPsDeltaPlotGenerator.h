#if !(defined VECPSDELTAPLOTGENERATOR)
#define VECPSDELTAPLOTGENERATOR

#include <vector>
#include <string>

#include "AMPTOOLS_AMPS/Vec_ps_refl.h"

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class VecPsDeltaPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum Hist_index{ kVecPsMass, kCosTheta, kPhi, 
                   kCosThetaH, kPhiH, kProd_Ang, kProdOffset, kt, 
                   kRecoilMass, kProtonPsMass, kRecoilPsMass, 
                   kLambda, kDalitz, kPhi_ProdVsPhi, kPhiOffsetVsPhi,
                   kCosThetaLV, kPhiLV, kPhi_ProdLVVsPhiLV,
                   kPhiOffsetLVVsPhiLV, kCosThetaPhi, kCosThetaHPhiH,
                   kCosThetaPhiLV, kTwoPsMass, kVecPsLVMass,
                   kVecPsPsMass, kPhiProdUVLV, kNumHists};

  VecPsDeltaPlotGenerator( const FitResults& results, Option opt);
  VecPsDeltaPlotGenerator( const FitResults& results );
  VecPsDeltaPlotGenerator( );

    static std::string numToString(Hist_index kHistName){
    switch(kHistName){
      case VecPsDeltaPlotGenerator::kVecPsMass: return "MVecPs"; 
      case VecPsDeltaPlotGenerator::kCosTheta: return "CosTheta";
      case VecPsDeltaPlotGenerator::kPhi: return "Phi";
      case VecPsDeltaPlotGenerator::kCosThetaH: return "CosTheta_H";
      case VecPsDeltaPlotGenerator::kPhiH: return "Phi_H";
      case VecPsDeltaPlotGenerator::kProd_Ang: return "Prod_Ang";
      case VecPsDeltaPlotGenerator::kProdOffset: return "ProdOffset";
      case VecPsDeltaPlotGenerator::kt: return "t";
      case VecPsDeltaPlotGenerator::kRecoilMass: return "MRecoil";
      case VecPsDeltaPlotGenerator::kProtonPsMass: return "MProtonPs";
      case VecPsDeltaPlotGenerator::kRecoilPsMass: return "MRecoilPs";
      case VecPsDeltaPlotGenerator::kLambda: return "Lambda";
      case VecPsDeltaPlotGenerator::kDalitz: return "Dalitz";
      case VecPsDeltaPlotGenerator::kPhi_ProdVsPhi: return "Phi_ProdVsPhi";
      case VecPsDeltaPlotGenerator::kPhiOffsetVsPhi: return "PhiOffsetVsPhi";
      case VecPsDeltaPlotGenerator::kCosThetaLV: return "CosThetaLV";
      case VecPsDeltaPlotGenerator::kPhiLV: return "PhiLV";
      case VecPsDeltaPlotGenerator::kPhi_ProdLVVsPhiLV: return "Phi_ProdLVVsPhiLV";
      case VecPsDeltaPlotGenerator::kPhiOffsetLVVsPhiLV: return "PhiOffsetLVVsPhiLV";
      case VecPsDeltaPlotGenerator::kCosThetaPhi: return "CosThetaPhi";
      case VecPsDeltaPlotGenerator::kCosThetaHPhiH: return "CosThetaHPhiH";
      case VecPsDeltaPlotGenerator::kCosThetaPhiLV: return "CosThetaPhiLV";
      case VecPsDeltaPlotGenerator::kTwoPsMass: return "MTwoPs";
      case VecPsDeltaPlotGenerator::kVecPsLVMass: return "MVecPsLV";
      case VecPsDeltaPlotGenerator::kVecPsPsMass: return "MVecPsPs";
      case VecPsDeltaPlotGenerator::kPhiProdUVLV: return "PhiProdUVLV";
      // Add more variables here if needed
            default: return "Unknown";
    }
}
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );

  std::map<std::string, Vec_ps_refl::VecPsReflArgs> m_args;
  
  void cacheArgs();
 
};

#endif

