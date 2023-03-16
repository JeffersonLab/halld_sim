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
	GDouble bigPhi = GPU_UVARS(4);
	GDouble phi = GPU_UVARS(5);

	GDouble W = 3.0*(0.5 - rho011)*sinSqTheta + rho011*(1.0 + 3.0*cosTheta*cosTheta) - 2.0*sqrt(3.0)*rho031*cos(phi)*sin2Theta - 2.0*sqrt(3.0)*rho03m1*cos(2.0*phi)*sinSqTheta; 

	W -= Pgamma*cos(2.0*bigPhi)*(3.0*rho133*sinSqTheta + rho111*(1.0 + 3.0*cosTheta*cosTheta) - 2.0*sqrt(3.0)*rho131*cos(phi)*sin2Theta - 2.0*sqrt(3.0)*rho13m1*cos(2.0*phi)*sinSqTheta);

	W -= Pgamma*sin(2.0*bigPhi)*(2.0*sqrt(3.0)*rho231*sin(phi)*sin2Theta + 2.0*sqrt(3.0)*rho23m1*sin(2.0*phi)*sinSqTheta);

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

