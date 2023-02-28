
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

////////////////////////////////////////////////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUVec_ps_refl_kernel( GPU_AMP_PROTO, int m_j, int m_m, int m_l, int m_r, int m_s, GDouble polAngle, GDouble polFraction )
{
	int iEvent = GPU_THIS_EVENT;

  	GDouble cosTheta =  GPU_UVARS(0);
	GDouble Phi =  GPU_UVARS(1);
	GDouble cosThetaH = GPU_UVARS(2);
	GDouble PhiH = GPU_UVARS(3);
	GDouble prod_angle = GPU_UVARS(4);
	GDouble MX = GPU_UVARS(5);
	GDouble MVec = GPU_UVARS(6);
	GDouble MPs = GPU_UVARS(7);
	
	///////////////////////////////////////////////////////////////////////////////////////////

	WCUComplex amplitude = CZero;

	for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
		//CPU --> clebschGordan(j1,j2,m1,m2,J,M) || GPU --> clebsch(j1,m1,j2,m2,J,M);
		GDouble hel_amp = clebsch(m_l, 0, 1, lambda, m_j, lambda);
		amplitude += Conjugate(wignerD( m_j, m_m, lambda, cosTheta, Phi )) * hel_amp * Conjugate(wignerD( 1, lambda, 0, cosThetaH, PhiH ));
  	} 
  
	GDouble Factor = sqrt(1 + m_s * polFraction);
	WCUComplex zjm = CZero;
	WCUComplex rotateY = { G_COS(  -1. * (prod_angle + polAngle*DegToRad) ) , G_SIN( -1. * (prod_angle + polAngle*DegToRad) ) }; // prod_angle and polAngle should be added

	if (m_r == 1)
		zjm = (amplitude * rotateY).m_dRe;
	if (m_r == -1) 
		zjm = (amplitude * rotateY).m_dIm;
		
  	Factor *= barrierFactor(MX, m_l, MVec, MPs);	

  	pcDevAmp[iEvent] = zjm * Factor;
}

void
GPUVec_ps_refl_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_j, int m_m, int m_l, int m_r, int m_s, GDouble polAngle, GDouble polFraction )

{  

  	GPUVec_ps_refl_kernel<<< dimGrid, dimBlock >>>
    		( GPU_AMP_ARGS, m_j, m_m, m_l, m_r, m_s, polAngle, polFraction);

}

