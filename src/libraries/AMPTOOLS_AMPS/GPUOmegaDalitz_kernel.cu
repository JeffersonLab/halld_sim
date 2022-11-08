
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUOmegaDalitz_kernel( GPU_AMP_PROTO, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta )
{
	int iEvent = GPU_THIS_EVENT;

	GDouble dalitz_z = GPU_UVARS(0);
	GDouble dalitz_sin3theta = GPU_UVARS(1);
	GDouble dalitz_phi = GPU_UVARS(2);
	
	///////////////////////////////////////////////////////////////////////////////////////////

	// dalitz parameters for 3-body vector decay
	GDouble G = G_SQRT( G_FABS(dalitz_phi * (1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * G_POW(dalitz_z,3/2.) * dalitz_sin3theta + 2 * dalitz_gamma * G_POW(dalitz_z,2) + 2 * dalitz_delta * G_POW(dalitz_z,5/2.) * dalitz_sin3theta)) );

  	pcDevAmp[iEvent] = G;
}

void
GPUOmegaDalitz_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta )

{  

  	GPUOmegaDalitz_kernel<<< dimGrid, dimBlock >>>
    		( GPU_AMP_ARGS, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta);

}

