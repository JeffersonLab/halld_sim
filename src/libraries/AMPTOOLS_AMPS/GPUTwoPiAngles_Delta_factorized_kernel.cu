
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
//
GPUTwoPiAngles_Delta_factorized_kernel( GPU_AMP_PROTO,
GDouble rho000,
GDouble rho100,
GDouble rho1m10,
GDouble rho111,
GDouble rho101,
GDouble rho1m11,
GDouble rho102,
GDouble rho1m12,
GDouble delta_rho011,
GDouble delta_rho031,
GDouble delta_rho03m1,
GDouble rho001
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
  GDouble cosTheta_piplus     = GPU_UVARS(5);
  GDouble sinSqTheta_piplus   = GPU_UVARS(6);
  GDouble sin2Theta_piplus    = GPU_UVARS(7);
  GDouble phi_piplus          = GPU_UVARS(8);
  GDouble bigPhi              = GPU_UVARS(9);
  // flag for pol (0.0 = no pol, 1.0 = pol)
  
  // GDouble rho001_val;
  // rho001_val = 2*(delta_rho111 + delta_rho133 - rho111); 
// Use beam assy. to eliminate one parameter

  GDouble W0_meson = 3.*(0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta_pim*cosTheta_pim - sqrt(2.)*rho100*sin2Theta_pim*cos(phi_pim) - rho1m10*sinSqTheta_pim*cos(2.*phi_pim)); 
  GDouble W1_meson= 3.*(rho111*sinSqTheta_pim + rho001*cosTheta_pim*cosTheta_pim - sqrt(2.)*rho101*sin2Theta_pim*cos(phi_pim) - rho1m11*sinSqTheta_pim*cos(2.*phi_pim)); 
  GDouble W2_meson = 3. * (sqrt(2.)*rho102*sin2Theta_pim*sin(phi_pim) + rho1m12*sinSqTheta_pim*sin(2.*phi_pim)); 

  GDouble W0_baryon = 3.*(0.5 - delta_rho011)*sinSqTheta_piplus + delta_rho011*(1.+3.*cosTheta_piplus*cosTheta_piplus) - 2.*sqrt(3.)*delta_rho031*cos(phi_piplus)*sin2Theta_piplus - 2.*sqrt(3.)*delta_rho03m1*cos(2.*phi_piplus)*sinSqTheta_piplus;
  // GDouble W1_baryon =  (3.*delta_rho133*sinSqTheta_piplus + delta_rho111*(1.+3.*cosTheta_piplus*cosTheta_piplus) - 2.*sqrt(3.)*delta_rho131*cos(phi_piplus)*sin2Theta_piplus - 2.*sqrt(3.)*delta_rho13m1*cos(2.*phi_piplus)*sinSqTheta_piplus);
  // GDouble W2_baryon = (2.*sqrt(3.)*delta_rho231*sin(phi_piplus)*sin2Theta_piplus + 2.*sqrt(3.)*delta_rho23m1*sin(2.*phi_piplus)*sinSqTheta_piplus);

  // GDouble Wpol = W0_meson*W0_baryon - Pgamma * cos(2.0 * bigPhi) * W1_meson * W1_baryon - Pgamma * sin(2.0 * bigPhi) * W2_meson * W2_baryon;

  /// Just another tests!!!!!! --- IGNORE ---

  // GDouble Wpol = W0_meson*W0_baryon -( Pgamma * cos(2.0 * bigPhi) * W1_meson + Pgamma * sin(2.0 * bigPhi) * W2_meson) * (W1_baryon + W2_baryon);
//Just another test
  GDouble Wpol = W0_meson*W0_baryon - Pgamma * cos(2.0 * bigPhi) *( W1_meson*W0_baryon) - Pgamma * sin(2.0 * bigPhi) *( W2_meson*W0_baryon);
  
  // GDouble Wpol = (W0_meson - Pgamma * cos(2.0 * bigPhi) * W1_meson - Pgamma * sin(2.0 * bigPhi) * W2_meson )*(W0_baryon - Pgamma * cos(2.0 * bigPhi) * W1_baryon - Pgamma * sin(2.0 * bigPhi) * W2_baryon);

  Wpol *= 1/(4.*PI);
  
  WCUComplex amp = { sqrt(fabs(Wpol)), 0.0 };
  pcDevAmp[iEvent] = amp;

}


void
GPUTwoPiAngles_Delta_factorized_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
  GDouble rho000,
  GDouble rho100,
  GDouble rho1m10,
  GDouble rho111,
  GDouble rho101,
  GDouble rho1m11,
  GDouble rho102,
  GDouble rho1m12,
  GDouble delta_rho011,
  GDouble delta_rho031,
  GDouble delta_rho03m1,
  GDouble rho001
)
{  

  GPUTwoPiAngles_Delta_factorized_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS,
    rho000,
    rho100,
    rho1m10,
    rho111,
    rho101,
    rho1m11,
    rho102,
    rho1m12,
    delta_rho011,
    delta_rho031,
    delta_rho03m1,
    rho001
    );
}
