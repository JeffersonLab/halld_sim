#include <cmath>
#include <complex>
#include "breakupMomentumComplex.h"

// mass0 = mass of parent
// mass1 = mass of first daughter
// mass2 = mass of second daughter

complex <GDouble> breakupMomentumComplex( GDouble mass0, GDouble mass1, GDouble mass2 ){
	
  complex <GDouble> q = std::sqrt(    mass0*mass0*mass0*mass0 + 
						  mass1*mass1*mass1*mass1 +
						  mass2*mass2*mass2*mass2 -
						  2.0*mass0*mass0*mass1*mass1 -
						  2.0*mass0*mass0*mass2*mass2 -
						  2.0*mass1*mass1*mass2*mass2  ) / (2.0 * mass0);

   return q;
	
}
