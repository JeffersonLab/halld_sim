
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUPhaseOffset_kernel(GPU_AMP_PROTO, GDouble phase)
{

  WCUComplex ans = {G_COS(phase), G_SIN(phase)};  
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUPhaseOffset_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble m_phase)
{
  GPUPhaseOffset_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, m_phase);
}
