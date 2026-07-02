#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUDeltaAngles_Rho_kernel( GPU_AMP_PROTO, GDouble delta_rho011, GDouble delta_rho031, GDouble delta_rho03m1,
                        GDouble delta_rho111, GDouble delta_rho133, GDouble delta_rho131, GDouble delta_rho13m1,
                        GDouble delta_rho231, GDouble delta_rho23m1, GDouble polAngle )
{
	int iEvent = GPU_THIS_EVENT;

	GDouble Pgamma = GPU_UVARS(0);
	GDouble cosTheta = GPU_UVARS(1);
	GDouble sinSqTheta = GPU_UVARS(2);
	GDouble sin2Theta = GPU_UVARS(3);
	GDouble bigPhi = GPU_UVARS(4);
	GDouble phi = GPU_UVARS(5);

	GDouble W = 3.0*(0.5 - delta_rho011)*sinSqTheta + delta_rho011*(1.0 + 3.0*cosTheta*cosTheta) - 2.0*sqrt(3.0)*delta_rho031*cos(phi)*sin2Theta - 2.0*sqrt(3.0)*delta_rho03m1*cos(2.0*phi)*sinSqTheta; 

	W -= Pgamma*cos(2.0*bigPhi)*(3.0*delta_rho133*sinSqTheta + delta_rho111*(1.0 + 3.0*cosTheta*cosTheta) - 2.0*sqrt(3.0)*delta_rho131*cos(phi)*sin2Theta - 2.0*sqrt(3.0)*delta_rho13m1*cos(2.0*phi)*sinSqTheta);

	W -= Pgamma*sin(2.0*bigPhi)*(2.0*sqrt(3.0)*delta_rho231*sin(phi)*sin2Theta + 2.0*sqrt(3.0)*delta_rho23m1*sin(2.0*phi)*sinSqTheta);

	W *= 1.0/(4.0*PI);

	WCUComplex amp = { sqrt( fabs( W ) ), 0 };

	pcDevAmp[iEvent] = amp;

}

void
GPUDeltaAngles_Rho_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                        GDouble delta_rho011, GDouble delta_rho031, GDouble delta_rho03m1,
                        GDouble delta_rho111, GDouble delta_rho133, GDouble delta_rho131, GDouble delta_rho13m1,
                        GDouble delta_rho231, GDouble delta_rho23m1, GDouble polAngle )

{

	GPUDeltaAngles_Rho_kernel<<< dimGrid, dimBlock >>>
		( GPU_AMP_ARGS, delta_rho011, delta_rho031, delta_rho03m1, delta_rho111, delta_rho133, delta_rho131, delta_rho13m1, delta_rho231, delta_rho23m1, polAngle );

}