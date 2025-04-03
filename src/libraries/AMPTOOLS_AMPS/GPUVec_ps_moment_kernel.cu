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


    // be careful to index using the proper integer corresponding to the enumeration in 
    // the C++ header file
    GDouble beamPolFraction = GPU_UVARS(0);
    GDouble cosTheta = GPU_UVARS(1);
    GDouble phi = GPU_UVARS(2);
    GDouble cosThetaH = GPU_UVARS(3);
    GDouble phiH = GPU_UVARS(4);
    GDouble prodAngle = GPU_UVARS(5);

    GDouble cos2prodAngle = G_COS(2*prodAngle);
    GDouble sin2prodAngle = G_SIN(2*prodAngle);
    GDouble theta = G_ACOS(cosTheta) * 180./PI;
    GDouble thetaH = G_ACOS(cosThetaH) * 180./PI;

    GDouble total = 0.0; // initialize the total "amplitude" to zero
    for(int i = 0; i < numberOfMoments; i++) {
      moment mom = moments[i];
      int Galpha = mom.alpha;
      int GJv = mom.Jv;
      int GLambda = mom.Lambda;
      int GJ = mom.J;
      int GM = mom.M;
      
      GDouble angle = (2*GJ + 1) / (4*PI ) * (2*GJv + 1) / (4*PI );
      if(Galpha == 0) {
        angle *= (
          wignerDSmall( GJ, GM, GLambda, theta ) * 
          wignerDSmall( GJv, GLambda, 0, thetaH ) * 
          G_COS( GM*phi + GLambda*phiH )
        );
      }
      else if(Galpha == 1) {
        angle *= -1.0 * beamPolFraction * cos2prodAngle * 
          wignerDSmall( GJ, GM, GLambda, theta ) * 
          wignerDSmall( GJv, GLambda, 0, thetaH ) * 
          G_COS( GM*phi + GLambda*phiH );
      }
      else if(Galpha == 2) {
        angle *= -1.0 * beamPolFraction * sin2prodAngle * 
          wignerDSmall( GJ, GM, GLambda, theta ) * 
          wignerDSmall( GJv, GLambda, 0, thetaH ) * 
          G_SIN( GM*phi + GLambda*phiH );
      }

      total += angle * H[i];
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