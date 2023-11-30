
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"
 
__global__ void
GPUTwoPSMoment_kernel(GPU_AMP_PROTO, GDouble *H, int *alpha, int *L, int *M, int nMoments )
{
    int iEvent = GPU_THIS_EVENT;

    GDouble pGamma = GPU_UVARS(0);
    GDouble cosTheta = GPU_UVARS(1);
    GDouble phi = GPU_UVARS(2);
    GDouble bigPhi = GPU_UVARS(3);

    GDouble total = 0;
    for(int imom = 0; imom < nMoments; imom++) { 
	   	   
	    int Galpha = alpha[imom];
	    int GL = L[imom];
	    int GM = M[imom];
	    
	    GDouble mom = 2.0 * sqrt( (2*GL + 1) / (4*PI ) );
	    if(Galpha == 0)
		   mom *= Y( GL, GM, cosTheta, phi ).m_dRe;
	    else if(Galpha == 1) 
		    mom *= pGamma * cos(2*bigPhi) * Y( GL, GM, cosTheta, phi ).m_dRe;
	    else if(Galpha == 2) 
		    mom *= -1 * pGamma * sin(2*bigPhi) * Y( GL, GM, cosTheta, phi ).m_dIm;
		    
	    // m = 0 only non-zero for alpha = 0, 1 but half the size of other m-projections
	    if(GM == 0 && Galpha < 2) mom *= 0.5;
		    
	    total += H[imom]*mom;
    }   

    WCUComplex amp = { sqrt( fabs( total ) ), 0 };

    pcDevAmp[iEvent] = amp;
}

void
GPUTwoPSMoment_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, int* alpha, int* L, int *M, int nMoments)
{

  // allocate memory and pass moment parameter array to GPU
  GDouble* d_H;
  int *d_alpha, *d_L, *d_M;
  cudaMalloc((void**)&d_H, nMoments * sizeof(GDouble));
  cudaMalloc((void**)&d_alpha, nMoments * sizeof(int));
  cudaMalloc((void**)&d_L, nMoments * sizeof(int));
  cudaMalloc((void**)&d_M, nMoments * sizeof(int));
  cudaMemcpy(d_H, &H[0], nMoments * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_alpha, &alpha[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_L, &L[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_M, &M[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );

  GPUTwoPSMoment_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_H, d_alpha, d_L, d_M, nMoments);

  cudaDeviceSynchronize();
  cudaFree(d_H);
}
