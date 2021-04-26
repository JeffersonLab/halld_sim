#include <stdio.h>
#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "breakupMomentum.cuh"
#include "barrierFactor.cuh"

#include "cuda.h"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

//////////////////////////////////////////// Useful Definitions ////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUVec_ps_refl_kernel( GPU_AMP_PROTO, int m_j, int m_m, int m_l, int m_r, int m_s, int m_3pi,
            GDouble polAngle, GDouble polFraction, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta)
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

///////////////////////////////////////////////////////////////////////////////////////////

   WCUComplex amplitude = CZero;

 // dalitz parameters for 3-body vector decay
  GDouble G = 1; // not relevant for 2-body vector decays
  if(m_3pi) G = G_SQRT(1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * G_POW(dalitz_z,3/2.) * dalitz_sin3theta + 2 * dalitz_gamma * G_POW(dalitz_z,2) + 2 * dalitz_delta * G_POW(dalitz_z,5/2.) * dalitz_sin3theta );

//  complex <GDouble> amplitude(0,0);
//  complex <GDouble> i(0,1);

  for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
	  //CPU --> clebschGordan(j1,j2,m1,m2,J,M) || GPU --> clebsch(j1,m1,j2,m2,J,M);
	  GDouble hel_amp = clebsch(m_l, 0, 1, lambda, m_j, lambda);
	  amplitude += Conjugate(wignerD( m_j, m_m, lambda, cosTheta, Phi )) * hel_amp * Conjugate(wignerD( 1, lambda, 0, cosThetaH, PhiH )) * G;
  } 
  
  GDouble Factor = sqrt(1 + m_s * polfrac);
  WCUComplex zjlambda = CZero;
  WCUComplex rotateY = { G_COS(  -1. * prod_angle ) , G_SIN( -1. * prod_angle ) };

  if (m_r == 1)
	  zjlambda = (amplitude * rotateY).m_dRe;
  if (m_r == -1) 
	  zjlambda = ic * (amplitude * rotateY).m_dIm;

  pcDevAmp[iEvent] = Factor * zjlambda;

}

void
GPUVec_ps_refl_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
	int m_j, int m_m, int m_l, int m_r, int m_s, int m_3pi, GDouble polAngle, GDouble polFraction,
         GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta)

{  

  GPUVec_ps_refl_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, m_j, m_m, m_l, m_r, m_s, m_3pi, polAngle, polFraction,
         dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta);

}
