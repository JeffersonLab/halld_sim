
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

////////////////////////////////////////////////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUSinglePS_kernel( GPU_AMP_PROTO, int m_r, int m_s )
{
	int iEvent = GPU_THIS_EVENT;

	GDouble prod_angle = GPU_UVARS(0);
	GDouble polFraction = GPU_UVARS(1);
	GDouble polAngle = GPU_UVARS(2);
	
	///////////////////////////////////////////////////////////////////////////////////////////

	WCUComplex amplitude = CZero;

  
	GDouble Factor = sqrt(1 + m_s * polFraction);
	WCUComplex rotateY = { G_COS(  -1. * (prod_angle + polAngle*DegToRad) ) , G_SIN( -1. * (prod_angle + polAngle*DegToRad) ) }; // prod_angle and polAngle should be added

	if (m_r == 1)
		amplitude = ( rotateY ).m_dRe;
	if (m_r == -1) 
		amplitude = ic * ( rotateY ).m_dIm;
		
  	pcDevAmp[iEvent] = amplitude * Factor;
}

void
GPUSinglePS_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_r, int m_s )

{  

  	GPUSinglePS_kernel<<< dimGrid, dimBlock >>>
    		( GPU_AMP_ARGS, m_r, m_s );

}

