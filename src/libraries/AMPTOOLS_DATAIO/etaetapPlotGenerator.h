#if !(defined ETAETAPPLOTGENERATOR)
#define ETAETAPPLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

//class FitResults;
class Kinematics;

// inheriting stuff from the PlotGenerator class here
class etaetapPlotGenerator : public PlotGenerator
{

public:

  enum Hist_index{
    khm12 = 0, khm13, khm23, khm1, khm2, khm3, kdltz, cosT, phiAng, PhiT, cosT_m23, Omega, cosT_phi, cosT_Phi, cosT_lab, phiAng_lab, cosT_m23_lab, phi_m23_lab,
    kNumHists
  }; //this is various histogram pointers
  
  etaetapPlotGenerator( const FitResults& results, Option opt);
  etaetapPlotGenerator( const FitResults& results ); //this is where we define the argument of our executable which is .fit file 
  etaetapPlotGenerator( );


  static std::string numToString(Hist_index kHistName){
    switch(kHistName){
      //adding one histogram name here for the testing
      case etaetapPlotGenerator::khm23: return "MPsPs";

      // case VecPsPlotGenerator::kVecPsMass: return "MVecPs"; 
      // case VecPsPlotGenerator::kCosTheta: return "CosTheta";
      // case VecPsPlotGenerator::kPhi: return "Phi";
      // case VecPsPlotGenerator::kCosThetaH: return "CosTheta_H";
      // case VecPsPlotGenerator::kPhiH: return "Phi_H";
      // case VecPsPlotGenerator::kProd_Ang: return "Prod_Ang";
      // case VecPsPlotGenerator::kt: return "t";
      // case VecPsPlotGenerator::kRecoilPsMass: return "MRecoilPs";
      // case VecPsPlotGenerator::kLambda: return "Lambda";
      // case VecPsPlotGenerator::kProdOffset: return "ProdOffset";
      // case VecPsPlotGenerator::kPhi_ProdVsPhi: return "Phi_ProdVsPhi";
      // case VecPsPlotGenerator::kPhiOffsetVsPhi: return "PhiOffsetVsPhi";
      // Add more variables here if needed
      default: return "Unknown";
    }
  }

private:

  double polAngle; // polarization angle in DEGREES
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );

  void createHistograms( );

};

#endif
