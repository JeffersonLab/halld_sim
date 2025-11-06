#if !(defined TWOPSPLOTGENERATOR)
#define TWOPSPLOTGENERATOR

#include <map>
#include <string>
#include <vector>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

// Forward declarations
class FitResults;
class Kinematics;

class TwoPsPlotGenerator : public PlotGenerator {
public:
  // Enumeration of histogram indices for the different observables
  enum HistogramIndex {
    k2PsMass = 0,
    kLambdaKMass,
    kLambdaPiMass,
    kLambdaMass,
    kdaughter1Mass, // Proton mass
    kdaughter2Mass, // Pi- mass
    kPiCosTheta,
    kPhiK,
    kPhiPi,
    kPhiLambda,
    kThetaK,
    kThetaPi,
    kThetaLambda,
    kMomK,
    kMomPi,
    kMomLambda,
    kPhi_LAB,
    kPhi,
    kphi,
    kPsi,
    kt,
    kCosThetaX_Lambda,
    kCosThetaY_Lambda,
    kCosThetaZ_Lambda,
    kCosThetaX_LambdaHel,
    kCosThetaY_LambdaHel,
    kCosThetaZ_LambdaHel,
    kCosTheta_LambdaHel,
    kPhi_LambdaHel,
    kphi_LambdaHel,
    kPsi_LambdaHel,
    kNumHists // Total number of histograms
  };

  // Constructors
  explicit TwoPsPlotGenerator(const FitResults &results);
  TwoPsPlotGenerator();

  // Projects a kinematic event onto histograms
  void projectEvent(Kinematics *kin);
  void projectEvent(Kinematics *kin, const string &reactionName);

  virtual ~TwoPsPlotGenerator();

private:
  // Method to initialize histograms
  void createHistograms();

  // Map to store reaction-specific angle data
  map<string, double> m_reactionAngleMap;
};

#endif // TWOPSPLOTGENERATOR
