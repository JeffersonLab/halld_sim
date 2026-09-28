#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUDeltaAngles_kernel( GPU_AMP_PROTO, GDouble rho011, GDouble rho031, GDouble rho03m1,
                        GDouble rho111, GDouble rho133, GDouble rho131, GDouble rho13m1,
                        GDouble rho231, GDouble rho23m1, GDouble polAngle )
{
	int iEvent = GPU_THIS_EVENT;

	GDouble Pgamma = GPU_UVARS(0);
	GDouble cosTheta = GPU_UVARS(1);
	GDouble sinSqTheta = GPU_UVARS(2);
	GDouble sin2Theta = GPU_UVARS(3);
  GDouble cosPhi = GPU_UVARS(4);
  GDouble cos2Phi = GPU_UVARS(5);
  GDouble sinPhi = GPU_UVARS(6);
  GDouble sin2Phi = GPU_UVARS(7);
  GDouble cos2BigPhi = GPU_UVARS(8);
  GDouble sin2BigPhi = GPU_UVARS(9);
  GDouble sqrt3 = sqrt(3.0);

  GDouble W = 3.*(0.5 - rho011)*sinSqTheta + rho011*(1.+3.*cosTheta*cosTheta) - 2.*sqrt3*rho031*cosPhi*sin2Theta - 2.*sqrt3*rho03m1*cos2Phi*sinSqTheta;
  
  W -= Pgamma*cos2BigPhi * (3.*rho133*sinSqTheta + rho111*(1.+3.*cosTheta*cosTheta) - 2.*sqrt3*rho131*cosPhi*sin2Theta - 2.*sqrt3*rho13m1*cos2Phi*sinSqTheta);
  
  W -= Pgamma*sin2BigPhi * (2.*sqrt3*rho231*sinPhi*sin2Theta + 2.*sqrt3*rho23m1*sin2Phi*sinSqTheta);

	W *= 1.0/(4.0*PI);

	WCUComplex amp = { sqrt( fabs( W ) ), 0 };

	pcDevAmp[iEvent] = amp;

}

void
GPUDeltaAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                        GDouble rho011, GDouble rho031, GDouble rho03m1,
                        GDouble rho111, GDouble rho133, GDouble rho131, GDouble rho13m1,
                        GDouble rho231, GDouble rho23m1, GDouble polAngle )

{

	GPUDeltaAngles_kernel<<< dimGrid, dimBlock >>>
		( GPU_AMP_ARGS, rho011, rho031, rho03m1, rho111, rho133, rho131, rho13m1, rho231, rho23m1, polAngle );

}

