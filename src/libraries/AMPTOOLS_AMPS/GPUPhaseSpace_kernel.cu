
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUPhaseSpace_kernel(GPU_AMP_PROTO, GDouble phase)
{

  WCUComplex ans = {G_Cos(phase), G_Sin(phase)};  
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUPhaseSpace_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble m_phase)
{
  GPUPhaseSpace_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, m_phase);
}
