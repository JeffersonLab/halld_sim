
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUComplexCoeff_kernel(GPU_AMP_PROTO, GDouble real, GDouble imag)
{
  WCUComplex ans = { real, imag };
  pcDevAmp[GPU_THIS_EVENT] = ans;
}

void
GPUComplexCoeff_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble real, GDouble imag)
{
  GPUComplexCoeff_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, real, imag);
}
