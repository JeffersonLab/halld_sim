
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUComplexCoeff_kernel(GPU_AMP_PROTO, GDouble param1, GDouble param2, bool represReIm)
{

  if(represReIm) {
    WCUComplex ans = { param1, param2 };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }
  else {
    WCUComplex ans = { param1*cos(param2), param1*sin(param2) };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }

}

void
GPUComplexCoeff_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble m_param1, GDouble m_param2, bool m_represReIm)
{
  GPUComplexCoeff_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, m_param1, m_param2, m_represReIm);
}
