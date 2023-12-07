
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"
 
__global__ void
GPUTwoPSMoment_kernel(GPU_AMP_PROTO, GDouble *H, int *indices, int nMoments, int maxL)
{
    int iEvent = GPU_THIS_EVENT;

    GDouble pGamma = GPU_UVARS(0);
    GDouble cosTheta = GPU_UVARS(1);
    GDouble phi = GPU_UVARS(2);
    GDouble bigPhi = GPU_UVARS(3);

    GDouble cos2bigPhi = G_COS(2*bigPhi);
    GDouble sin2bigPhi = G_SIN(2*bigPhi);
    GDouble theta = G_ACOS(cosTheta) * 180.0 / PI;

    // compute required wignerDSmall values
    GDouble wigner[100];
    for(int iL = 0; iL <= maxL; iL++) {
	for(int iM = 0; iM <= iL; iM++) {
	    wigner[iL*10 + iM] = wignerDSmall( iL, iM, 0, theta );
	}
    }

    GDouble total = 0;
    for(int imom = 0; imom < nMoments; imom++) { 

	    int Galpha = indices[imom] / 100;
	    int GL = indices[imom] / 10 % 10;
	    int GM = indices[imom] % 10;
	    int GLM = indices[imom] % 100;

	    GDouble mom = 2.0 * (2*GL + 1) / (4*PI);
	
	    // compute moments via wignerDSmall
	    if(Galpha == 0)
	    	   mom *= wigner[GLM] * cos( GM*phi ); // Y( GL, GM, cosTheta, phi ).m_dRe;  
	    else if(Galpha == 1) 
		   mom *= pGamma * cos2bigPhi * wigner[GLM] * cos( GM*phi );
	    else if(Galpha == 2) 
		   mom *= -1.0 * pGamma * sin2bigPhi * wigner[GLM] * sin( GM*phi );
	    
	    // m = 0 only non-zero for alpha = 0, 1 but half the size of other m-projections
	    if(GM == 0 && Galpha < 2) mom *= 0.5;
		    
	    total += H[imom]*mom;
    }   

    WCUComplex amp = { G_SQRT( G_FABS( total ) ), 0 };

    pcDevAmp[iEvent] = amp;

}

void
GPUTwoPSMoment_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble* H, int *indices, int nMoments, int maxL)
{

  // allocate memory and pass moment parameter array to GPU
  GDouble* d_H;
  int *d_indices;
  cudaMalloc((void**)&d_H, nMoments * sizeof(GDouble));
  cudaMalloc((void**)&d_indices, nMoments * sizeof(int));
  cudaMemcpy(d_H, &H[0], nMoments * sizeof(GDouble), cudaMemcpyHostToDevice );
  cudaMemcpy(d_indices, &indices[0], nMoments * sizeof(int), cudaMemcpyHostToDevice );

  GPUTwoPSMoment_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS, d_H, d_indices, nMoments, maxL);

  cudaDeviceSynchronize();
  cudaFree(d_H);
  cudaFree(d_indices);
}
