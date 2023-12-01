
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"
 
__global__ void
GPUTwoPSMoment_kernel(GPU_AMP_PROTO, GDouble *H, int *alpha, int *S, int *Lambda, int *J, int *M, int nMoments )
{
    int iEvent = GPU_THIS_EVENT;

    GDouble pGamma = GPU_UVARS(0);
    GDouble cosTheta = GPU_UVARS(1);
    GDouble phi = GPU_UVARS(2);
    GDouble cosThetaH = GPU_UVARS(1);
    GDouble phiH = GPU_UVARS(2);
    GDouble bigPhi = GPU_UVARS(3);
    GDouble theta = acos(cosTheta) * 180./PI;
    GDouble thetaH = acos(cosThetaH) * 180./PI;

    GDouble total = 0.0;
    for(int imom = 0; imom < nMoments; imom++) { 
	   	   
	    int Galpha = alpha[imom];
	    int GS = S[imom];
	    int GLambda = Lambda[imom];
	    int GJ = J[imom];
	    int GM = M[imom];
	    
	    GDouble mom = 2.0 * sqrt( (2*GJ + 1) / (4*PI ) ) * sqrt( (2*GS + 1) / (4*PI ) );
	    if(Galpha == 0)
		   mom *= wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
	    else if(Galpha == 1) 
		    mom *= -1.0 * pGamma * cos(2*bigPhi) * wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
	    else if(Galpha == 2) 
		    mom *= -1.0 * pGamma * sin(2*bigPhi) * wignerDSmall( GJ, GM, GLambda, theta ) * wignerDSmall( GS, GLambda, 0, thetaH ) * cos( GM*phi + GLambda*phiH );
		    
	    total += H[imom]*mom;
    }   

    //GDouble total = totalAmp.m_dRe;

    WCUComplex amp = { sqrt( fabs( total ) ), 0 };

    pcDevAmp[iEvent] = amp;
}

void
GPUTwoPSMoment_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, int* alpha, int *S, int *Lambda, int *J, int *M, int nMoments)
{

  // allocate memory and pass moment parameter array to GPU
  GDouble* d_H;
  int *d_alpha, *d_S, *d_Lambda, *d_J, *d_M;
  cudaMalloc((void**)&d_H, nMoments * sizeof(GDouble));
  cudaMalloc((void**)&d_alpha, nMoments * sizeof(int));
  cudaMalloc((void**)&d_S, nMoments * sizeof(int));
  cudaMalloc((void**)&d_Lambda, nMoments * sizeof(int));
  cudaMalloc((void**)&d_J, nMoments * sizeof(int));
  cudaMalloc((void**)&d_M, nMoments * sizeof(int));

  cudaMemcpy(d_H, &H[0], nMoments * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_alpha, &alpha[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_S, &S[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_Lambda, &Lambda[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_J, &J[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );
  cudaMemcpy(d_M, &M[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );

  GPUTwoPSMoment_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_H, d_alpha, d_S, d_Lambda, d_J, d_M, nMoments);

  cudaDeviceSynchronize();
  cudaFree(d_H);
}
