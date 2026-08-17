#if !(defined TWOPSPLOTGENERATOR)
#define TWOPSPLOTGENERATOR

#include <map>
#include <string>
#include <vector>

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

class Kinematics;

class TwoPsPlotGenerator : public PlotGenerator
{

public:

  // Particle labels (TLatex) injected here, indexed to kin->particle(i):
  //   [0] = beam, [1] = Lambda, [2] = Ks, [3] = pi+
  TwoPsPlotGenerator( const FitResults& results, Option opt,
                      const std::vector<std::string>& partName
                       = { "beam", "1", "2", "3" } );

  enum Hist_index{
    // 'kFirsthist' ties the first entry in the enum to the plotter (in this case TwoPsPlotter.cc).
    kFirstHist = 0,
    // invariant masses / Dalitz
    kMass12 = 0, kMass13, kMass23, kMass1, kMass2, kMass3, kDalitz,
    // lab frame
    kCosTheta_lab, kPhi_lab, kCosTheta_m23_lab, kPhi_m23_lab,
    // helicity frame -- polarization angle (event-global) + zero-angle diagnostic
    kBigPhi_hel, kBigPhiOffset_hel,
    // helicity frame -- Mass 2 (Ks)
    kCosThetaM2_hel, kPhiM2_hel, kPhiMinusBigPhiM2_hel, kCosThetaM2_m23_hel, kCosThetaM2_Phi_hel,
    // helicity frame -- Mass 3 (pi+)
    kCosThetaM3_hel, kPhiM3_hel, kPhiMinusBigPhiM3_hel, kCosThetaM3_m23_hel, kCosThetaM3_Phi_hel,
    // helicity frame -- phi vs Phi correlation (real + offset diagnostic)
    kPhiM2_BigPhi_hel, kPhiM2_BigPhiOffset_hel,
    kNumHists
  };

  static std::string numToString( Hist_index i );

private:

  // Cached per-reaction polarization info (filled once by cacheArgs()).
  struct PolInfo {
    bool   polInPhotonP4 = false;   // 5-arg Zlm: pol stored in beam photon px/py
    double polAngleDeg   = 0.0;     // 7-arg Zlm: fixed polarization angle (degrees)
  };

  std::map<std::string, PolInfo> m_polInfo;   // key = reaction name
  std::vector<std::string>       m_partName;  // display labels, indexed as above

  void cacheArgs();
  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const std::string& reactionName );

};

#endif
