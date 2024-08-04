
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
	
  	GDouble alpha_term = GPU_UVARS(0);
  	GDouble beta_term = GPU_UVARS(1);
  	GDouble gamma_term = GPU_UVARS(2);
  	GDouble delta_term = GPU_UVARS(3);
	GDouble lambda = GPU_UVARS(4);
	
	// dalitz parameters for 3-body vector decay

	GDouble G = G_SQRT( G_FABS( lambda *  ( 1 + dalitz_alpha * alpha_term + dalitz_beta * beta_term + dalitz_gamma * gamma_term + dalitz_delta * delta_term ) ) );

  	pcDevAmp[iEvent] = G;
}

void
GPUOmegaDalitz_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta )

{  

  	GPUOmegaDalitz_kernel<<< dimGrid, dimBlock >>>
    		( GPU_AMP_ARGS, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta);

}

