
#if !defined(DECAYANGLES)
#define DECAYANGLES

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

using std::complex;
using namespace std;

vector< TVector3 > getHelicityAxes(const TVector3 &parentCM, const TVector3 &inverseCM);

vector< TVector3 > getGJAxes(const TVector3 &parentCM, const TVector3 &inverseCM, const TVector3 &inverseParent);

double getPhiProd(double polAngle, const TLorentzVector &parentLab, const TLorentzVector &beamLab, const TLorentzVector &targetLab, int whichFrame, bool upperVertex);
double getPhiProd(const TVector3 &eps, const TLorentzVector &parentLab, const TLorentzVector &beamLab, const TLorentzVector &targetLab, int whichFrame, bool upperVertex);

vector< double > getOneStepAngles(const TLorentzVector &parentLab, const TLorentzVector &daughterLab, const TLorentzVector &beamLab, const TLorentzVector &targetLab, int whichFrame, bool upperVertex = true);

vector< double > getTwoStepAngles(const TLorentzVector &parentLab, const TLorentzVector &daughterLab, const TLorentzVector &granddaughter1Lab, const TLorentzVector &granddaughter2Lab, const TLorentzVector &beamLab, const TLorentzVector &targetLab, int whichFrame, bool upperVertex);

#endif
