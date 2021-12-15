
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
 
__global__ void
GPUPiecewise_kernel(GPU_AMP_PROTO, GDouble * params1, GDouble * params2, int nBins, bool represReIm )
{

  int iEvent = GPU_THIS_EVENT;
  int tempBin = GPU_UVARS(0);

  // some thread debugging info
  //unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
  //if(threadIdx.x == 0){
  //  printf("Hello from block %d dim %d, thread %d: index = %d\n %d %f %f\n", blockIdx.x, blockDim.x, threadIdx.x, index, tempBin, params1[tempBin], params2[tempBin]);
  //}

  if(represReIm) {
    WCUComplex ans = { params1[tempBin], params2[tempBin] };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }
  else {
    WCUComplex ans = { params1[tempBin]*cos(params2[tempBin]), params1[tempBin]*sin(params2[tempBin]) };
    pcDevAmp[GPU_THIS_EVENT] = ans;
  }
}

void
GPUPiecewise_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* params1, GDouble* params2, int nBins, bool represReIm)
{

  // allocate memory and pass piecewise parameter array to GPU
  double* d_params1;
  double* d_params2;
  cudaMalloc((void**)&d_params1, nBins * sizeof(double));
  cudaMalloc((void**)&d_params2, nBins * sizeof(double));
  cudaMemcpy(d_params1, &params1[0], nBins * sizeof(double), cudaMemcpyHostToDevice );
  cudaMemcpy(d_params2, &params2[0], nBins * sizeof(double), cudaMemcpyHostToDevice );

  GPUPiecewise_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_params1, d_params2, nBins, represReIm);
}
