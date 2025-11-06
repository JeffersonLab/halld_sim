#include <math.h>
#include <cuda.h>
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

__global__ void
// GPUKStarHyperon_kernel( GPU_AMP_PROTO,
//   GDouble alpha, GDouble Sigma, GDouble Ox, GDouble P,
//   GDouble T, GDouble Oz, GDouble polAngle )
GPUKStarHyperon_kernel( GPU_AMP_PROTO,
                    GDouble alpha, GDouble rho111, GDouble rho001, GDouble Ox, GDouble P,
                    GDouble T, GDouble Oz, GDouble polAngle ) {

  int iEvent = GPU_THIS_EVENT;

  GDouble Pgamma    = GPU_UVARS(0);
  GDouble cosThetaX = GPU_UVARS(1);
  GDouble cosThetaY = GPU_UVARS(2);
  GDouble cosThetaZ = GPU_UVARS(3);
  GDouble phi       = polAngle * 0.017453293 + GPU_UVARS(4);

  GDouble I = 1.0 + alpha * cosThetaY * P;
  //I -= Pgamma * cos(2.0 * phi) * (Sigma + alpha * cosThetaY * T);
  I -= Pgamma * cos(2.0 * phi) * (2 * rho111 + rho001 + alpha * cosThetaY * T);
  I += Pgamma * sin(2.0 * phi) * alpha * (cosThetaX * Ox + cosThetaZ * Oz);

  WCUComplex amp = { sqrt(fabs(I)), 0 };
  pcDevAmp[iEvent] = amp;
}


// void
// GPUKStarHyperon_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
//                   GDouble alpha, GDouble Sigma, GDouble Ox, GDouble P,
//                   GDouble T, GDouble Oz, GDouble polAngle ) {

//   GPUKStarHyperon_kernel<<< dimGrid, dimBlock >>>(
//     GPU_AMP_ARGS, alpha, Sigma, Ox, P, T, Oz, polAngle );
// }
void
GPUKStarHyperon_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                  GDouble alpha, GDouble rho111, GDouble rho001, GDouble Ox, GDouble P,
                  GDouble T, GDouble Oz, GDouble polAngle ) {

  GPUKStarHyperon_kernel<<< dimGrid, dimBlock >>>(
    GPU_AMP_ARGS, alpha, rho111, rho001, Ox, P, T, Oz, polAngle );
}