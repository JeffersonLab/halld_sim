
#if !defined(HELICITYFRAME)
#define HELICITYFRAME

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

#include "GPUManager/GPUCustomTypes.h"

using std::complex;
using namespace std;

vector <double> getthetaphi(TLorentzVector particle1, TLorentzVector z, TLorentzVector InverseOfX, TLorentzVector zrf, TLorentzVector particle2, TLorentzVector xrf);

vector <double> getthetaphi(TLorentzVector particle1, TLorentzVector z, TLorentzVector InverseOfX, TLorentzVector zrf);

#endif
