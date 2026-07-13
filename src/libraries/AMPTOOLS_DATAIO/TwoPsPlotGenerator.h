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
    // ---------- SDME A ----------
    kA000, kA100, kA1m10, kA111, kA001, kA101, kA1m11, kA102, kA1m12,
    // ---------- Polarization B ----------
    kB1, kB2, kB3, kB4, kB5,
    // ---------- Extended Ax ----------
    kPx0, kAx000, kAx10c0, kAx10s0, kAx1m1c0, kAx1m1s0,
    kPx1, kAx001, kAx10c1, kAx10s1, kAx1m1c1, kAx1m1s1,
    kPx2, kAx002, kAx10c2, kAx10s2, kAx1m1c2, kAx1m1s2,
    // ---------- Extended Ay ----------
    kPy0, kAy000, kAy10c0, kAy10s0, kAy1m1c0, kAy1m1s0,
    kPy1, kAy001, kAy10c1, kAy10s1, kAy1m1c1, kAy1m1s1,
    kPy2, kAy002, kAy10c2, kAy10s2, kAy1m1c2, kAy1m1s2,
    // ---------- Extended Az ----------
    kPz0, kAz000, kAz10c0, kAz10s0, kAz1m1c0, kAz1m1s0,
    kPz1, kAz001, kAz10c1, kAz10s1, kAz1m1c1, kAz1m1s1,
    kPz2, kAz002, kAz10c2, kAz10s2, kAz1m1c2, kAz1m1s2,
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
