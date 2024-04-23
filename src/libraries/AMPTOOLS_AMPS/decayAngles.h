
#if !defined(DECAYANGLES)
#define DECAYANGLES

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

using std::complex;
using namespace std;

vector< TVector3 > getHelicityAxes(TVector3 parentCM, TVector3 inverseCM);

vector< TVector3 > getGJAxes(TVector3 parentCM, TVector3 inverseCM, TVector3 inverseParent);

double getPhiProd(double polAngle, TLorentzVector parentLab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex);

vector< double > getOneStepAngles(TLorentzVector parentLab, TLorentzVector daughterLab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex = true);

vector< double > getTwoStepAngles(TLorentzVector parentLab, TLorentzVector daughterLab, TLorentzVector granddaughter1Lab, TLorentzVector granddaughter2Lab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex);

#endif
