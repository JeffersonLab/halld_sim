
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
#include "DblReggeHelper_FastPi.cuh"
//#include "AMPTOOLS_AMPS/DblReggeHelper.cuh"

__global__ void
DblRegge_FastPi_kernel(GPU_AMP_PROTO, GDouble S0, GDouble b_pi, int charge){

	int iEvent = GPU_THIS_EVENT;

	// here we need to be careful to index the user-defined
	// data with the proper integer corresponding to the
	// enumeration in the C++ header file

	//user vars as defined in enum in header:


	GDouble s12 = GPU_UVARS(0);
	GDouble s23 = GPU_UVARS(1);
	GDouble t1 = GPU_UVARS(2);
//	GDouble t2 = GPU_UVARS(3);
	GDouble s = GPU_UVARS(4);
	GDouble u3 = GPU_UVARS(5);
	GDouble beamM2 = GPU_UVARS(6);
	GDouble p1M2 = GPU_UVARS(7);
	GDouble p2M2 = GPU_UVARS(8);
	GDouble recoilM2 = GPU_UVARS(9);
//	GDouble up1 = GPU_UVARS(10);
//	GDouble up2 = GPU_UVARS(11);



	WCUComplex amp =  GPU_calcAmplitude(s, s12, s23, t1, u3, S0, b_pi,beamM2, p1M2, p2M2, recoilM2, charge);

//if((amp.Re() + amp.Im()) > 35)
//{
//printf( "amp: %f \n", amp);
//printf( "u3: %f \n", u3);
//printf( "t1: %f \n", t1);
//printf( "s23: %f \n", s23);
//printf( "s12: %f \n", s12);
//}	

	pcDevAmp[iEvent] = amp;
}



void GPUDblRegge_FastPi_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		GDouble S0, GDouble b_pi, int charge  )
{

	DblRegge_FastPi_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS, S0, b_pi, charge );
}


