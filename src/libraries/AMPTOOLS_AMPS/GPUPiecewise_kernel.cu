
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUPiecewise_kernel(GPU_AMP_PROTO)
{

  int tempBin = GPU_UVARS(0);
  WCUComplex ans = { m_paramsRe[tempBin], m_paramsIm[tempBin] };
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUPiecewise_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO)
{
  GPUPiecewise_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS);
}
