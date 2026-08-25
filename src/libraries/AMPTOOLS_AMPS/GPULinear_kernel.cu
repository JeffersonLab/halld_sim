#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPULinear_kernel( GPU_AMP_PROTO, GDouble real_p0, GDouble real_p1,
                  GDouble imag_p0, GDouble imag_p1  ){
    
  int iEvent = GPU_THIS_EVENT;

  GDouble mass = GPU_UVARS(0);
  
  GDouble real_tot = real_p0 + real_p1 * mass;
  GDouble imag_tot = imag_p0 + imag_p1 * mass;
  
  WCUComplex amp = { real_tot, imag_tot };
  pcDevAmp[iEvent] = amp;
}

void
GPULinear_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                GDouble real_p0, GDouble real_p1,
                GDouble imag_p0, GDouble imag_p1 ){
  
  GPULinear_kernel<<< dimGrid, dimBlock >>>
   ( GPU_AMP_ARGS, real_p0, real_p1, imag_p0, imag_p1 );
}

