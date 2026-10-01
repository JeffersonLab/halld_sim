
#if !defined(DECAYANGLESTEST)
#define DECAYANGLESTEST

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

using std::complex;
using namespace std;

typedef struct {
    double bigPhiX;
    double cosThetaX;
    double phiX;
    double cosThetaH;
    double phiH;
    double lambda;
} VecPSAngles;

typedef struct {
    double bigPhiLV;
    double cosThetaLV;
    double phiLV;
} LVAngles;

void calcVectorDecayAngles( VecPSAngles& outputAngles,
                            const TLorentzVector& vecVPS,       // four-vectors in the vec-ps frame
                            const TLorentzVector& vecDaught1VPS,
                            const TLorentzVector& vecDaught2VPS,
                            const TVector3& zAxisParent );

LVAngles getLVAnglesGJ( const TLorentzVector& beamLab,
                        const TLorentzVector& recoilLab,
                        const TLorentzVector& protonLab );

LVAngles getLVAnglesHelicity( const TLorentzVector& beamLab,
                              const TLorentzVector& recoilLab,
                              const TLorentzVector& protonLab );

VecPSAngles getVecPSAnglesGJ( const TLorentzVector& beamLab,
                              const TLorentzVector& vecPSLab,
                              const TLorentzVector& vecLab,
                              const TLorentzVector& vecDaught1Lab,
                              const TLorentzVector& vecDaught2Lab );

VecPSAngles getVecPSAnglesHelicity( const TLorentzVector& beamLab,
                                    const TLorentzVector& vecPSLab,
                                    const TLorentzVector& vecLab,
                                    const TLorentzVector& vecDaught1Lab,
                                    const TLorentzVector& vecDaught2Lab );

#endif
