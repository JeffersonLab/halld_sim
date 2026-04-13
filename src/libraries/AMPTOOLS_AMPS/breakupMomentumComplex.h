#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"



#if !defined(BREAKUPMOMENTUMCOMPLEX)
#define BREAKUPMOMENTUMCOMPLEX

// mass0 = mass of parent
// mass1 = mass of first daughter
// mass2 = mass of second daughter

complex <GDouble> breakupMomentumComplex( GDouble mass0, GDouble mass1, GDouble mass2 );

#endif
