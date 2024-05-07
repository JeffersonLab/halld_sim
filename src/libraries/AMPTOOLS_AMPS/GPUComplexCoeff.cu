
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUComplexCoeff_kernel(GPU_AMP_PROTO, GDouble m_value_re, GDouble m_value_im)
{
  WCUComplex ans = { m_value_re, m_value_im };
  pcDevAmp[GPU_THIS_EVENT] = ans;
}

void
GPUComplexCoeff_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble m_value_re, GDouble m_value_im)
{
  GPUComplexCoeff_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, m_value_re, m_value_im);
}
