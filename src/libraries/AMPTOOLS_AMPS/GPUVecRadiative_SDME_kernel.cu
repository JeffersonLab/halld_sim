
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUVecRadiative_SDME_kernel( GPU_AMP_PROTO, GDouble rho000, GDouble rho100,
                          GDouble rho1m10, GDouble rho111, GDouble rho001,
                          GDouble rho101, GDouble rho1m11, GDouble rho102,
                          GDouble rho1m12 ){

  int iEvent = GPU_THIS_EVENT;

  // here we need to be careful to index the user-defined
  // data with the proper integer corresponding to the
  // enumeration in the C++ header file

  GDouble polFrac = GPU_UVARS(0);
  GDouble cosTheta = GPU_UVARS(1);
  GDouble sinSqTheta = GPU_UVARS(2);
  GDouble sin2Theta = GPU_UVARS(3);
  GDouble bigPhi = GPU_UVARS(4);
  GDouble phi = GPU_UVARS(5);


  GDouble rho110 = 0.5*(1-rho000);

  GDouble W = 1 - sinSqTheta*rho110 - cosTheta*cosTheta*rho000 +
              sinSqTheta*cos(2*phi)*rho1m10 +
              sqrt(2.)*rho100*sin2Theta*cos(phi);
	
  W -= polFrac * cos(2*bigPhi) * ( 2*rho111 + sinSqTheta*(rho001-rho111) +
                                   sinSqTheta*cos(2*phi)*rho1m11 +
                                   sqrt(2.)*rho101*sin2Theta*cos(phi) );
	
  W += polFrac * sin(2*bigPhi) * ( rho1m12*sinSqTheta*sin(2*phi) +
                                   sqrt(2.)*rho102*sin2Theta*sin(phi) );
	
  W *= 3/(8*PI);

  WCUComplex amp = { sqrt( fabs( W ) ), 0 };

  pcDevAmp[iEvent] = amp;
}


void
GPUVecRadiative_SDME_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                           GDouble rho000, GDouble rho100, GDouble rho1m10,
                           GDouble rho111, GDouble rho001, GDouble rho101,
                           GDouble rho1m11, GDouble rho102, GDouble rho1m12 )
{  

  GPUVecRadiative_SDME_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, rho000, rho100, rho1m10, rho111, rho001,
      rho101, rho1m11, rho102, rho1m12 );
}
