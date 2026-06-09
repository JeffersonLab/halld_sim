#if !defined(BREAKUPMOMENTUM)
#define BREAKUPMOMENTUM

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

// mass0 = mass of parent
// mass1 = mass of first daughter
// mass2 = mass of second daughter

double breakupMomentum( double mass0, double mass1, double mass2 );

std::complex<GDouble> breakupMomentumComplex( GDouble mass0, GDouble mass1, GDouble mass2 );

#endif
