
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
//
GPUTwoPiAngles_Delta_DoubleSDMEs_unpol_kernel( GPU_AMP_PROTO, 
  //order is like the appearing order in the paper
  //Wmm
  GDouble r00_33_0,  //GDouble r00_11_0,
  GDouble r11_33_0,  GDouble r11_11_0,
  GDouble r1m1_33_0, GDouble r1m1_11_0,
  GDouble r10_33_0,  GDouble r10_11_0,

  //Wbar
  GDouble r00_31_0,  GDouble r00_3m1_0,
  GDouble r11_31_0,  GDouble r11_3m1_0,
  GDouble r1m1_31_0, GDouble r1m1_3m1_0,
  GDouble r10_31_0,  GDouble r10_3m1_0,

  //Wtil
  GDouble rt1m1_31_0, GDouble rt1m1_3m1_0,
  GDouble rt10_31_0,  GDouble rt10_3m1_0
)
  {

  int iEvent = GPU_THIS_EVENT;

  // here we need to be careful to index the user-defined
  // data with the proper integer corresponding to the
  // enumeration in the C++ header file


  GDouble Pgamma              = GPU_UVARS(0);
  GDouble cosTheta_pim        = GPU_UVARS(1);
  GDouble sinSqTheta_pim      = GPU_UVARS(2);
  GDouble sin2Theta_pim       = GPU_UVARS(3);
  GDouble phi_pim             = GPU_UVARS(4);
  GDouble cosTheta_proton     = GPU_UVARS(5);
  GDouble sinSqTheta_proton   = GPU_UVARS(6);
  GDouble sin2Theta_proton    = GPU_UVARS(7);
  GDouble phi_proton          = GPU_UVARS(8);
  GDouble bigPhi              = GPU_UVARS(9);
  // flag for pol (0.0 = no pol, 1.0 = pol)
  
  GDouble sinSqTh = sinSqTheta_proton;
  GDouble sin2Th = sin2Theta_proton;
  GDouble cosSqTh = 1.0 - sinSqTh;

  GDouble cphi = cos(phi_proton);
  GDouble sphi = sin(phi_proton);
  GDouble c2phi = cos(2.0 * phi_proton);
  GDouble s2phi = sin(2.0 * phi_proton);

  GDouble sinSqThPi = sinSqTheta_pim;
  GDouble sin2ThPi = sin2Theta_pim;
  GDouble cosSqThPi = 1.0 - sinSqThPi;

  GDouble cphiPi = cos(phi_pim);
  GDouble sphiPi = sin(phi_pim);
  GDouble c2phiPi = cos(2.0 * phi_pim);
  GDouble s2phiPi = sin(2.0 * phi_pim);

  GDouble r00_11_0 = 0.5 - r00_33_0 - 2.0 * r11_11_0 - 2.0 * r11_33_0;

  GDouble W00_0 = r00_33_0 * sinSqTh + r00_11_0 * (1.0/3.0 + cosSqTh);
  GDouble Wb00_0 = (2.0 / sqrt(3.0)) * (r00_31_0 * cphi * sin2Th + r00_3m1_0 * c2phi * sinSqTh);

  GDouble W11_0 = r11_33_0 * sinSqTh + r11_11_0 * (1.0/3.0 + cosSqTh);
  GDouble Wb11_0 = (2.0 / sqrt(3.0)) * (r11_31_0 * cphi * sin2Th + r11_3m1_0 * c2phi * sinSqTh);

  GDouble W1m1_0 = r1m1_33_0 * sinSqTh + r1m1_11_0 * (1.0/3.0 + cosSqTh);
  GDouble Wb1m1_0 = (2.0 / sqrt(3.0)) * (r1m1_31_0 * cphi * sin2Th + r1m1_3m1_0 * c2phi * sinSqTh);
  GDouble Wt1m1_0 = (2.0 / sqrt(3.0)) * (rt1m1_31_0 * sphi * sin2Th + rt1m1_3m1_0 * s2phi * sinSqTh);

  GDouble W10_0 = r10_33_0 * sinSqTh + r10_11_0 * (1.0/3.0 + cosSqTh);
  GDouble Wb10_0 = (2.0 / sqrt(3.0)) * (r10_31_0 * cphi * sin2Th + r10_3m1_0 * c2phi * sinSqTh);
  GDouble Wt10_0 = (2.0 / sqrt(3.0)) * (rt10_31_0 * sphi * sin2Th + rt10_3m1_0 * s2phi * sinSqTh);
  GDouble W0 = cosSqThPi * (W00_0 - Wb00_0) + sinSqThPi * (W11_0 - Wb11_0) - sinSqThPi * (c2phiPi * (W1m1_0 - Wb1m1_0) + s2phiPi * Wt1m1_0) - sqrt(2.0) * sin2ThPi * (cphiPi * (W10_0 - Wb10_0) + sphiPi * Wt10_0);

  WCUComplex amp = { sqrt(fabs(W0)), 0.0 };
  pcDevAmp[iEvent] = amp;

}


void
GPUTwoPiAngles_Delta_DoubleSDMEs_unpol_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
  GDouble r00_33_0,  //GDouble r00_11_0,
  GDouble r11_33_0,  GDouble r11_11_0,
  GDouble r1m1_33_0, GDouble r1m1_11_0,
  GDouble r10_33_0,  GDouble r10_11_0,
  GDouble r00_31_0,  GDouble r00_3m1_0,
  GDouble r11_31_0,  GDouble r11_3m1_0,
  GDouble r1m1_31_0, GDouble r1m1_3m1_0,
  GDouble r10_31_0,  GDouble r10_3m1_0,
  GDouble rt1m1_31_0, GDouble rt1m1_3m1_0,
  GDouble rt10_31_0,  GDouble rt10_3m1_0
)
{  

  GPUTwoPiAngles_Delta_DoubleSDMEs_unpol_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, 
      r00_33_0, //r00_11_0, 
      r11_33_0, r11_11_0, 
      r1m1_33_0, r1m1_11_0, 
      r10_33_0, r10_11_0,
      r00_31_0, r00_3m1_0, 
      r11_31_0, r11_3m1_0, 
      r1m1_31_0, r1m1_3m1_0, 
      r10_31_0, r10_3m1_0,
      rt1m1_31_0, rt1m1_3m1_0, 
      rt10_31_0, rt10_3m1_0
    );
}
