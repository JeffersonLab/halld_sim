
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUTwoLeptonAngles_kernel( GPU_AMP_PROTO, GDouble rho000, GDouble rho100,
                           GDouble rho1m10, GDouble rho111, GDouble rho001,
		           GDouble rho101, GDouble rho1m11, GDouble rho102,
		           GDouble rho1m12, GDouble polAngle ){

  int iEvent = GPU_THIS_EVENT;

  // here we need to be careful to index the user-defined
  // data with the proper integer corresponding to the
  // enumeration in the C++ header file

  GDouble Pgamma = GPU_UVARS(0);
  GDouble cosTheta = GPU_UVARS(1);
  GDouble sinSqTheta = GPU_UVARS(2);
  GDouble sin2Theta = GPU_UVARS(3);
  GDouble bigPhi = polAngle*0.017453293 + GPU_UVARS(4);
  GDouble phi = GPU_UVARS(5);

  GDouble W = 0.5*(1. + rho000) - 0.5*(3.*rho000 - 1.)*cosTheta*cosTheta + sqrt(2.)*rho100*sin2Theta*cos(phi) + rho1m10*sinSqTheta*cos(2.*phi);
	
  W -= Pgamma*cos(2.*bigPhi) * (rho111*(1 + cosTheta*cosTheta) + rho001*sinSqTheta + sqrt(2.)*rho101*sin2Theta*cos(phi) + rho1m11*sinSqTheta*cos(2.*phi));
	
  W += Pgamma*sin(2.*bigPhi) * (sqrt(2.)*rho102*sin2Theta*sin(phi) + rho1m12*sinSqTheta*sin(2.*phi));
	
  W *= 3./(8.*PI);

  WCUComplex amp = { sqrt( fabs( W ) ), 0 };

  pcDevAmp[iEvent] = amp;
}


void
GPUTwoLeptonAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
                   GDouble rho000, GDouble rho100, GDouble rho1m10,
		   GDouble rho111, GDouble rho001, GDouble rho101,
		   GDouble rho1m11, GDouble rho102, GDouble rho1m12,
		   GDouble polAngle )
{  

  GPUTwoLeptonAngles_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, rho000, rho100, rho1m10, rho111, rho001,
      rho101, rho1m11, rho102, rho1m12, polAngle );
}
