#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

struct moment {
  char name[32]; // assuming max length of name is 32 characters
  int alpha;
  int Jv;
  int Lambda;
  int J;
  int M;
};

__global__ void
GPUVec_ps_moment_kernel(GPU_AMP_PROTO, GDouble *H, moment *moments, int numberOfMoments )
{
    int iEvent = GPU_THIS_EVENT;

    GDouble total = 0.0; // initialize the total "amplitude" to zero
    for(int i = 0; i < numberOfMoments; i++) {
      GDouble angular_distribution = GPU_UVARS(i);
      total += angular_distribution * H[i];
    }   

  // since AmpTools is hard coded to handle squared amplitudes, we return the square root of the total
  WCUComplex amp = { G_SQRT( G_FABS( total ) ), 0 };

  pcDevAmp[iEvent] = amp;
}

void
GPUVec_ps_moment_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, moment* moments, int numberOfMoments)
{
  // allocate memory and pass moment parameter array to GPU
  GDouble* d_H;
  moment *d_moments;
  cudaMalloc((void**)&d_H, numberOfMoments * sizeof(GDouble));
  cudaMalloc((void**)&d_moments, numberOfMoments * sizeof(moment));

  cudaMemcpy(d_H, &H[0], numberOfMoments * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_moments, &moments[0], numberOfMoments * sizeof(moment), cudaMemcpyHostToDevice );

  GPUVec_ps_moment_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_H, d_moments, numberOfMoments);

  cudaDeviceSynchronize();
  cudaFree(d_H);
  cudaFree(d_moments);
}