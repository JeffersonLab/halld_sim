#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

__global__ void GPUPiPiSWaveAMPK_kernel( GPU_AMP_PROTO )
{
	int iEvent = GPU_THIS_EVENT;

	WCUComplex ans = { GPU_UVARS(0), GPU_UVARS(1) };
    	pcDevAmp[iEvent] = ans;	
}

void GPUPiPiSWaveAMPK_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO )
{  

  	GPUPiPiSWaveAMPK_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS );

}