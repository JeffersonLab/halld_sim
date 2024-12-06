
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

__global__ void
GPUVec_ps_refl_kernel( GPU_AMP_PROTO )
{
	int iEvent = GPU_THIS_EVENT;

	WCUComplex ans = { GPU_UVARS(0), GPU_UVARS(1) };
    	pcDevAmp[iEvent] = ans;	
}

void
GPUVec_ps_refl_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO )

{  

  	GPUVec_ps_refl_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS );

}

