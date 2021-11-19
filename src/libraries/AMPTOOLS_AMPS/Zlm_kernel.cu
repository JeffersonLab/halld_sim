#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
#include "GPUUtils/wignerD.cuh"

__global__ void
Zlm_kernel( GPU_AMP_PROTO, int j, int m, int r, int s ){

  int iEvent = GPU_THIS_EVENT;

  // here we need to be careful to index the user-defined
  // data with the proper integer corresponding to the
  // enumeration in the C++ header file

  GDouble pGamma = GPU_UVARS(0);
  GDouble cosTheta = GPU_UVARS(1);
  GDouble phi = GPU_UVARS(2);
  GDouble bigPhi = GPU_UVARS(3);

  GDouble factor = sqrt(1 + s * pGamma);
  GDouble zlm = 0;

  WCUComplex rotateY = { G_COS( bigPhi ), -G_SIN( bigPhi ) };
  if( r == 1 ){
    zlm = ( Y( j, m, cosTheta, phi ) * rotateY ).Re();
  }
  if( r == -1 ){
    zlm = ( Y( j, m, cosTheta, phi ) * rotateY ).Im();
  }

  WCUComplex amp = { factor * zlm, 0 };

  pcDevAmp[iEvent] = amp;
}


void
GPUZlm_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
              int j, int m, int r, int s  )
{

  Zlm_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS, j, m, r, s );
}

