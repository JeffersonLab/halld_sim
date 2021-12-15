
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

__device__ double c_arrayRe[5];
__device__ double c_arrayIm[5];
 
__global__ void
GPUPiecewise_kernel(GPU_AMP_PROTO, GDouble * paramsRe, GDouble * paramsIm, int nBins)
{

  int iEvent = GPU_THIS_EVENT;
  int tempBin = GPU_UVARS(0);

  //extern __shared__ double bufRe[]; // size is not stated
  //extern __shared__ double bufIm[]; // size is not stated

  paramsRe[0] = c_arrayRe[0];
  paramsIm[0] = c_arrayIm[0];

  if(threadIdx.x == 100){
    printf("Hello from block %d dim %d, thread %d: index = %d \n", blockIdx.x, blockDim.x, threadIdx.x, blockIdx.x * blockDim.x + threadIdx.x);
    printf("%d %f\n", tempBin, paramsRe[0]);
  }
  //printf("%d %f\n", nBins, bufRe[0]);
  //printf("%f", paramsRe[blockIdx.x * blockDim.x + threadIdx.x]);
  //paramsIm[blockIdx.x * blockDim.x + threadIdx.x];

  //bufRe[threadIdx.x] = paramsRe[blockIdx.x * blockDim.x + threadIdx.x];
  //bufIm[nBins + threadIdx.x] = 1; //paramsIm[blockIdx.x * blockDim.x + threadIdx.x];

  WCUComplex ans = { 1, 0 }; //bufRe[tempBin], bufIm[tempBin] };
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUPiecewise_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* paramsRe, GDouble* paramsIm, int nBins)
{
  //cudaMemcpy(m_paramsReGPU, &paramsRe[0], paramsRe.size()*sizeof(GDouble), cudaMemcpyHostToDevice);
  //cudaMemcpy(m_paramsImGPU, &paramsIm[0], paramsIm.size()*sizeof(GDouble), cudaMemcpyHostToDevice);

  double* d_paramsRe = 0;
  double* d_paramsIm = 0;
  cudaMalloc((void**)&d_paramsRe, 5 * sizeof(double));
  cudaMalloc((void**)&d_paramsIm, 5 * sizeof(double));
  cudaMemcpyToSymbol("c_arrayRe", paramsRe, 5 * sizeof(double), 0, cudaMemcpyHostToDevice );
  cudaMemcpyToSymbol("c_arrayIm", paramsIm, 5 * sizeof(double), 0, cudaMemcpyHostToDevice );

  printf("%f %f %d %d %d\n", paramsRe[4], paramsIm[4], nBins, sizeof(double), 2*nBins*sizeof(double));
  GPUPiecewise_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_paramsRe, d_paramsIm, nBins);
}

/*
__global__ void kernel(float* d_array){ d_array[0] = c_array[0]; }
 
void test(){
    float* d_array = 0;
    float h_array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};   
    cudaMalloc((void**)&d_array, 10 * sizeof(float));
    cudaMemcpyToSymbol("c_array", h_array, sizeof(float)*10, 0, cudaMemcpyHostToDevice );
    kernel<<< 1, 1 >>>(d_array);
}
*/