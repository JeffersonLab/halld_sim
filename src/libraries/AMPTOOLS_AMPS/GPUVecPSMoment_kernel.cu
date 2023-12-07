
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"
 
__global__ void
GPUVecPSMoment_kernel(GPU_AMP_PROTO, GDouble *H, int *indices, int nMoments )
{
    int iEvent = GPU_THIS_EVENT;

    GDouble pGamma = GPU_UVARS(0);
    GDouble cosTheta = GPU_UVARS(1);
    GDouble phi = GPU_UVARS(2);
    GDouble cosThetaH = GPU_UVARS(1);
    GDouble phiH = GPU_UVARS(2);
    GDouble bigPhi = GPU_UVARS(3);

    GDouble cos2bigPhi = G_COS(2*bigPhi);
    GDouble sin2bigPhi = G_SIN(2*bigPhi);
    GDouble theta = G_ACOS(cosTheta) * 180./PI;
    GDouble thetaH = G_ACOS(cosThetaH) * 180./PI;

    GDouble total = 0.0;
    for(int imom = 0; imom < nMoments; imom++) { 

    	    int Galpha = indices[imom] / 10000;
	    int GS = indices[imom] / 1000 % 10;
	    int GLambda = indices[imom] / 100 % 10;
	    int GJ = indices[imom] / 10 % 10;
	    int GM = indices[imom] %10;
	    
	    GDouble mom = 2.0 * (2*GJ + 1) / (4*PI ) * (2*GS + 1) / (4*PI );
	    if(Galpha == 0)
		   mom *= wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
	    else if(Galpha == 1) 
		    mom *= -1.0 * pGamma * cos2bigPhi * wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
	    else if(Galpha == 2) 
		    mom *= -1.0 * pGamma * sin2bigPhi * wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
		    
	    total += H[imom]*mom;
    }   

    WCUComplex amp = { G_SQRT( G_FABS( total ) ), 0 };

    pcDevAmp[iEvent] = amp;
}

void
GPUVecPSMoment_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, int* indices, int nMoments)
{

  // allocate memory and pass moment parameter array to GPU
  GDouble* d_H;
  int *d_indices;
  cudaMalloc((void**)&d_H, nMoments * sizeof(GDouble));
  cudaMalloc((void**)&d_indices, nMoments * sizeof(int));

  cudaMemcpy(d_H, &H[0], nMoments * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_indices, &indices[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );

  GPUVecPSMoment_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_H, d_indices, nMoments);

  cudaDeviceSynchronize();
  cudaFree(d_H);
  cudaFree(d_indices);
}
