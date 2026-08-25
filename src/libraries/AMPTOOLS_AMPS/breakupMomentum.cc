#include <cmath>
#include <complex>
#include "breakupMomentum.h"

// mass0 = mass of parent
// mass1 = mass of first daughter
// mass2 = mass of second daughter

double breakupMomentum( double mass0, double mass1, double mass2 ){
	
        double q;

	// fabs -- correct?  consistent w/ previous E852 code
        q = sqrt( fabs(   mass0*mass0*mass0*mass0 + 
						  mass1*mass1*mass1*mass1 +
						  mass2*mass2*mass2*mass2 -
						  2.0*mass0*mass0*mass1*mass1 -
						  2.0*mass0*mass0*mass2*mass2 -
						  2.0*mass1*mass1*mass2*mass2  ) ) / (2.0 * mass0);

        return q;
	
}

std::complex<GDouble> breakupMomentumComplex( GDouble mass0, GDouble mass1, GDouble mass2 ){
	
  std::complex<GDouble> q = std::sqrt(    mass0*mass0*mass0*mass0 + 
						  mass1*mass1*mass1*mass1 +
						  mass2*mass2*mass2*mass2 -
						  2.0*mass0*mass0*mass1*mass1 -
						  2.0*mass0*mass0*mass2*mass2 -
						  2.0*mass1*mass1*mass2*mass2  ) / (2.0 * mass0);

   return q;
	
}
