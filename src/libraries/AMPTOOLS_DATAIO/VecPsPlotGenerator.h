#if !(defined VCPSPLOTGENERATOR)
#define VECPSPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class VecPsPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum Hist_index{ kVecPsMass, kCosTheta, kPhi, kCosThetaH, kPhiH, 
                   kProd_Ang, kProdOffset, kt, kRecoilMass, kProtonPsMass, 
                   kRecoilPsMass, kLambda, kDalitz, kPhi_ProdVsPhi, 
                   kPhiOffsetVsPhi, kNumHists};

  VecPsPlotGenerator( const FitResults& results, Option opt);
  VecPsPlotGenerator( const FitResults& results );
  VecPsPlotGenerator( );

    static std::string numToString(Hist_index kHistName){
    switch(kHistName){
      case VecPsPlotGenerator::kVecPsMass: return "MVecPs"; 
      case VecPsPlotGenerator::kCosTheta: return "CosTheta";
      case VecPsPlotGenerator::kPhi: return "Phi";
      case VecPsPlotGenerator::kCosThetaH: return "CosTheta_H";
      case VecPsPlotGenerator::kPhiH: return "Phi_H";
      case VecPsPlotGenerator::kProd_Ang: return "Prod_Ang";
      case VecPsPlotGenerator::kt: return "t";
      case VecPsPlotGenerator::kRecoilPsMass: return "MRecoilPs";
      case VecPsPlotGenerator::kLambda: return "Lambda";
      case VecPsPlotGenerator::kProdOffset: return "ProdOffset";
      case VecPsPlotGenerator::kPhi_ProdVsPhi: return "Phi_ProdVsPhi";
      case VecPsPlotGenerator::kPhiOffsetVsPhi: return "PhiOffsetVsPhi";
      // Add more variables here if needed
            default: return "Unknown";
    }
}
 
private:
  
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );
 
};

#endif

