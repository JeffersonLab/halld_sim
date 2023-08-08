
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
 
__global__ void
GPUPiecewise_kernel(GPU_AMP_PROTO, GDouble * params1, GDouble * params2, int nBins, bool represReIm )
{

  int iEvent = GPU_THIS_EVENT;
  
#ifdef AMPTOOLS_GDOUBLE_FP64
  long* tempBin = (long*)&(GPU_UVARS(0));
#else
  int* tempBin = (int*)&(GPU_UVARS(0));
 #endif

  if(represReIm) {
    WCUComplex ans = { params1[*tempBin], params2[*tempBin] };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }
  else {
    WCUComplex ans = { params1[*tempBin]*cos(params2[*tempBin]), params1[*tempBin]*sin(params2[*tempBin]) };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }
}

void
GPUPiecewise_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* params1, GDouble* params2, int nBins, bool represReIm)
{

  // allocate memory and pass piecewise parameter array to GPU
  GDouble* d_params1;
  GDouble* d_params2;
  cudaMalloc((void**)&d_params1, nBins * sizeof(GDouble));
  cudaMalloc((void**)&d_params2, nBins * sizeof(GDouble));
  cudaMemcpy(d_params1, &params1[0], nBins * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_params2, &params2[0], nBins * sizeof(GDouble), cudaMemcpyHostToDevice );

  GPUPiecewise_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_params1, d_params2, nBins, represReIm);

  cudaDeviceSynchronize();
  cudaFree(d_params1);
  cudaFree(d_params2);
}
