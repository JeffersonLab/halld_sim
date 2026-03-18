// The purpose of getVectorDecayAngles() is to return the helicity angles 
// of the normal to the decay plane of the vector parent particle.
// This function also returns the variable lambda used in the 
// omega->3pi decay amplitude.

// The purpose of getXDecayAngles() is to return the helicity angles of 
// particle X in the reaction gamma + p -> X + p.
// This function also returns the angle between the polarization
// vector and the normal to the production plane of particle X.
// Both functions return a vector of doubles containing the polar 
// angle theta, the azimuthal angle phi, and either lambda or bigPhi 
// The angles are calculated in the helicity frame of the parent 
// The code uses b1->omega pi and omega->3pi as an example,
// but can be used for any vector decay.

// For a two particle decay, the function arguments are the lab 
// frame TLorentzVectors of: the daughter particle, the parent 
// particle, the inverse of the X-axis and the reference Z-axis.

// For a three particle decay, the function arguments are the same 
// with the addition of a second daughter particle.
// Note that in this case the normal to the decay plane is calculated as
// the cross product between the first daughter and the second daughter vectors.

// The function could also return the decay angles in
// other frames when called with the appropriate argument change.
// Gottfried-Jackson (commented in code): The z-axis is equal to the
//        direction of the incoming beam in the particle X rest frame.
// Adair RF: The z-axis is equal to the direction of the incoming beam
//        in the center of mass system.
#if !defined(VECPSANGLES)
#define VECPSANGLES

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"

//#include "GPUManager/GPUCustomTypes.h"

using std::complex;
using namespace std;

constexpr double DEG_TO_RAD = TMath::Pi() / 180.0;

vector <double> getVectorDecayAngles( const TLorentzVector& beamPLab,
                                      const TLorentzVector& particleXLab,
                                      const TLorentzVector& parentLab, 
                                      const TLorentzVector& firstDaughterLab, 
                                      const TLorentzVector& secondDaughterLab);

vector<double> getXDecayAngles( double polAngle, 
                                const TLorentzVector& beamLab, 
                                const TLorentzVector& beamPLab,
                                const TLorentzVector& particleXLab,
                                const TLorentzVector& daughterLab);

#endif
