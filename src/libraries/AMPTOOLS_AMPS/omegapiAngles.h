
#if !defined(OMEGAPIANGLES)
#define OMEGAPIANGLES

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

//#include "GPUManager/GPUCustomTypes.h"

using std::complex;
using namespace std;

vector <double> getomegapiAngles(TLorentzVector daughter, TLorentzVector parent, TLorentzVector InverseOfX, TLorentzVector rf, TLorentzVector seconddaughter);

vector <double> getomegapiAngles(double polAngle, TLorentzVector daughter, TLorentzVector parent, TLorentzVector InverseOfX, TLorentzVector rf);

#endif
