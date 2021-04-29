
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

// Macro to ease definition of loops
#define LOOP(INDEX,START,END,INC) for (int INDEX=START;INDEX<=END;INDEX+=INC)
////////////////////////////////////////////////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUVec_ps_refl_kernel( GPU_AMP_PROTO, int m_j, int m_m, int m_l, int m_r, int m_s, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction )
{
	int iEvent = GPU_THIS_EVENT;

  	GDouble cosTheta =  GPU_UVARS(0);
	GDouble Phi =  GPU_UVARS(1);
	GDouble cosThetaH = GPU_UVARS(2);
	GDouble PhiH = GPU_UVARS(3);
	GDouble prod_angle = GPU_UVARS(4);
	GDouble polfrac = GPU_UVARS(5);
	GDouble dalitz_z = GPU_UVARS(6);
	GDouble dalitz_sin3theta = GPU_UVARS(7);
	GDouble M3Pi = GPU_UVARS(8);
	GDouble M4Pi = GPU_UVARS(9);

	///////////////////////////////////////////////////////////////////////////////////////////

	WCUComplex amplitude = CZero;

	GDouble G = G_SQRT( 1 + (2 * dalitz_alpha * dalitz_z) + (2 * dalitz_beta * G_POW(dalitz_z, (double) 3/2) * dalitz_sin3theta) + (2 * dalitz_gamma * G_POW(dalitz_z, (double) 2.0)) + (2 * dalitz_delta * G_POW(dalitz_z, (double) 5/2) * dalitz_sin3theta) );

	for (int lambda = -1; lambda <= 1; lambda++) { //omega helicity
		GDouble hel_amp = clebsch(m_l, 0, 1, lambda, m_j, lambda);
		//CPU --> clebschGordan(j1,j2,m1,m2,J,M) || GPU --> clebsch(j1,m1,j2,m2,J,M);
		
		amplitude += Conjugate(wignerD( m_j, m_m, lambda, cosTheta, Phi )) * hel_amp * Conjugate(wignerD( 1, lambda, 0, cosThetaH, PhiH )) * G;
	}//loop over lambda
		
	GDouble Factor = G_SQRT(1 + m_s * polfrac);
	WCUComplex zjm = CZero;
	WCUComplex rotateY = { G_COS(prod_angle), -1.*G_SIN(prod_angle) };
	WCUComplex product = amplitude * rotateY;
	if (m_r == 1)
        	zjm = product.Re();
	if (m_r == -1)
       		zjm = product.Im();

	// E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
	GDouble kinFactor = barrierFactor(M4Pi, m_l, M3Pi, 0.139);
  	//kinFactor *= sqrt(3.) * sqrt(2.*m_l + 1.);
  	Factor *= kinFactor;	

  	pcDevAmp[iEvent] = zjm * Factor;
}

void
GPUVec_ps_refl_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_j, int m_m, int m_l, int m_r, int m_s, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction)

{  

  	GPUVec_ps_refl_kernel<<< dimGrid, dimBlock >>>
    		( GPU_AMP_ARGS, m_j, m_m, m_l, m_r, m_s, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction);

}

