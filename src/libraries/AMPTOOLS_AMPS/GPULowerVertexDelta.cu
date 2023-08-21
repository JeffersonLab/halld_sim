#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

<<<<<<< HEAD
#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"
=======
#include "GPUUtils/wignerD.cuh"
>>>>>>> master
///////////////////////////////////////////////////////////////////////////////
__device__ WCUComplex CZero = { 0, 0 };
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPULowerVertexDelta_kernel( GPU_AMP_PROTO, int m_d, int m_p, int m_c, int m_s )
{
	int iEvent = GPU_THIS_EVENT;

	GDouble cosTheta = GPU_UVARS(0);
	GDouble phi = GPU_UVARS(1);
	
	GDouble lambda_Delta = m_d / 2.;
	GDouble lambda_proton = m_p / 2.;
	
	WCUComplex amp = CZero;

	if ( m_c == 1 )
		amp = 2/3.* m_s * wignerD( 3/2, lambda_Delta, lambda_proton, cosTheta, phi );
	if ( m_c == -1 )
		amp = 2/3.* m_s * Conjugate( wignerD( 3/2, lambda_Delta, lambda_proton, cosTheta, phi ) );

	pcDevAmp[iEvent] = amp;

}

void
GPULowerVertexDelta_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, int m_d, int m_p, int m_c, int m_s )

{

	GPULowerVertexDelta_kernel<<< dimGrid, dimBlock >>>
		( GPU_AMP_ARGS, m_d, m_p, m_c, m_s );

}

