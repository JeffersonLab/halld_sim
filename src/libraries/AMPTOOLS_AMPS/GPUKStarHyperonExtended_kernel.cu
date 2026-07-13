#include <math.h>
#include <cuda.h>
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

// GPU-safe POD (no std::string)
struct term_ijk {
  int i;
  int j;
  int k;
};

// __global__ void
// GPUKStarHyperonExtended_kernel(
//   GPU_AMP_PROTO,
//   int model,
//   GDouble alpha,
//   GDouble* coeffs,          // device array of coefficients (size nTerms)
//   term_ijk* terms,          // device array of (i,j,k) (size nTerms)
//   int nTerms,
//   bool schillingIncluded,
//   GDouble* rhos,            // device array size 9 if schillingIncluded else nullptr
//   GDouble polAngle
// ) {

//   int iEvent = GPU_THIS_EVENT;

//   GDouble Pgamma = GPU_UVARS(0);                          // kPgamma
//   GDouble bigPhi = polAngle * 0.017453293 + GPU_UVARS(1); // kBigPhi rotated

//   GDouble g_func[3];
//   g_func[0] = 1.0;
//   g_func[1] = -Pgamma * G_COS(2.0 * bigPhi);
//   g_func[2] = -Pgamma * G_SIN(2.0 * bigPhi);

//   // ----- extended cross-term sum -----
//   GDouble I_ext = 0.0;
//   if (model >= 2) {
//     for (int ii = 0; ii < 3; ++ii) {
//       const GDouble phi = coeffs[18 + ii] * 0.017453293;
//       const GDouble cph = G_COS(phi);
//       const GDouble sph = G_SIN(phi);

//       const GDouble g = g_func[ii];

//       const int jL = model;  
//       const int jT0 = (jL == 4) ? 2 : jL + 1;
//       const int jT1 = (jL == 2) ? 4 : jL - 1;

//       for (int kk = 0; kk < 6; ++kk) {
//         const int idxT = ii * 6 + kk;        // 0..17
//         const int idxL = 21 + ii * 6 + kk;   // 21..38

//         const GDouble T = coeffs[idxT];
//         const GDouble L = coeffs[idxL];

//         const GDouble b_k = GPU_UVARS(5 + kk);

//         // nx term uses T*cos(phi), ny uses T*sin(phi), nz uses L
//         I_ext += (T * cph) * g * GPU_UVARS(jT0) * b_k;
//         I_ext += (T * sph) * g * GPU_UVARS(jT1) * b_k;
//         I_ext += (L)       * g * GPU_UVARS(jL) * b_k;
//       }
//     }

//   } else {
//     for (int t = 0; t < nTerms; ++t) {
//       int ii = terms[t].i;
//       int jj = 2 + terms[t].j;
//       int kk = 5 + terms[t].k;

//       GDouble n_j = GPU_UVARS(jj);  // kN0..kN2
//       GDouble b_k = GPU_UVARS(kk);  // kB0..kB5
      
//       // if (threadIdx.x == 0 && blockIdx.x == 0) {
//       //   printf(
//       //       "t=%d ii=%d jj=%d kk=%d\n g_i=%g n_j=%g b_k=%g coeffs=%g\n "
//       //       "kN0=%g kN1=%g kN2=%g\n "
//       //       "kB0=%g kB1=%g kB2=%g kB3=%g kB4=%g kB5=%g\n",
//       //       t, ii, jj, kk,
//       //       g_func[ii], n_j, b_k, coeffs[t],
//       //       GPU_UVARS(2), GPU_UVARS(3), GPU_UVARS(4),
//       //       GPU_UVARS(5), GPU_UVARS(6), GPU_UVARS(7),
//       //       GPU_UVARS(8), GPU_UVARS(9), GPU_UVARS(10)
//       //   );
//       // }
//       I_ext += coeffs[t] * g_func[ii] * n_j * b_k;
//     }
//   }
//   I_ext *= 1.0 / (4.0 * PI);

//   // ----- optional Schilling–Wolf SDME -----
//   GDouble W = 1.0;
//   if (schillingIncluded) {
//     GDouble rho000  = rhos[0];
//     GDouble rho100  = rhos[1];
//     GDouble rho1m10 = rhos[2];
//     GDouble rho111  = rhos[3];
//     GDouble rho001  = rhos[4];
//     GDouble rho101  = rhos[5];
//     GDouble rho1m11 = rhos[6];
//     GDouble rho102  = rhos[7];
//     GDouble rho1m12 = rhos[8];

//     // basis shortcuts
//     GDouble cosSqTheta = GPU_UVARS(6); 
//     GDouble sinSqTheta = 1.0 - cosSqTheta;
//     GDouble b10c   = GPU_UVARS(7);     
//     GDouble b10s   = GPU_UVARS(8);     
//     GDouble b1m1c  = GPU_UVARS(9);     
//     GDouble b1m1s  = GPU_UVARS(10);    

//     const GDouble sqrt2 = G_SQRT(2.0);

//     W = 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosSqTheta - sqrt2*rho100*b10c - rho1m10*b1m1c;
//     W += g_func[1] * (rho111*sinSqTheta + rho001*cosSqTheta - sqrt2*rho101*b10c - rho1m11*b1m1c);
//     W += g_func[2] * (sqrt2*rho102*b10s + rho1m12*b1m1s);
//     W *= 3.0 / (4.0 * PI);
//   }

//   GDouble I_total =  W + alpha * I_ext;


//   WCUComplex amp = { G_SQRT(G_FABS(I_total)), 0 };
//   pcDevAmp[iEvent] = amp;
// }

__global__ void
GPUKStarHyperonExtended_kernel(
  GPU_AMP_PROTO,
  int model,
  GDouble alpha,
  GDouble* coeffs,          // device array of coefficients (size nTerms)
  term_ijk* terms,          // device array of (i,j,k) (size nTerms)
  int nTerms,
  bool schillingIncluded,
  GDouble* rhos,            // device array size 9 if schillingIncluded (size 17 for model 2) else nullptr
  GDouble polAngle
) {

  int iEvent = GPU_THIS_EVENT;

  GDouble Pgamma = GPU_UVARS(0);                          // kPgamma
  GDouble bigPhi = polAngle * 0.017453293 + GPU_UVARS(1); // kBigPhi rotated

  GDouble g_func[3];
  g_func[0] = 1.0;
  g_func[1] = -Pgamma * G_COS(2.0 * bigPhi);
  g_func[2] = -Pgamma * G_SIN(2.0 * bigPhi);
  // g_func[1] =  Pgamma * G_SIN(2.0 * bigPhi);
  // g_func[2] =  -Pgamma * G_COS(2.0 * bigPhi);
  // g_func[1] = 1.0;
  // g_func[2] = 1.0;

  // ----- extended cross-term sum -----
  GDouble H[3] = {};
  GDouble FH[3][2] = {
    {1.0, 0.0}, 
    {1.0, 0.0},
    {1.0, 0.0}
  }; // for factorized model


  if (model >= 3) {
    for (int ii = 0; ii < 3; ++ii) {
      const GDouble phi = coeffs[18 + ii] * 0.017453293;
      const GDouble cph = G_COS(phi);
      const GDouble sph = G_SIN(phi);

      const GDouble g = g_func[ii];

      const int jL = model - 1;  
      const int jT0 = (jL == 4) ? 2 : jL + 1;
      const int jT1 = (jL == 2) ? 4 : jL - 1;

      for (int kk = 0; kk < 6; ++kk) {
        const int idxT = ii * 6 + kk;        // 0..17
        const int idxL = 21 + ii * 6 + kk;   // 21..38

        const GDouble T = coeffs[idxT];
        const GDouble L = coeffs[idxL];

        const GDouble b_k = GPU_UVARS(5 + kk);

        // nx term uses T*cos(phi), ny uses T*sin(phi), nz uses L
        H[ii] += (T * cph) * g * GPU_UVARS(jT0) * b_k;
        H[ii] += (T * sph) * g * GPU_UVARS(jT1) * b_k;
        H[ii] += (L)       * g * GPU_UVARS(jL)  * b_k;
      }
    }

  } else if (model == 2) {
    for (int t = 0; t < nTerms; ++t) {
      int ii = terms[t].i;
      int jj = 2 + terms[t].j;
      int kk = 5 + terms[t].k;

      GDouble n_j = GPU_UVARS(jj);  // kN0..kN2
      GDouble b_k = GPU_UVARS(kk);  // kB0..kB5

      int signature = (terms[t].j == 1 ? 0 : 1);  // y -> 0, x/z -> 1
      FH[ii][signature] += alpha * coeffs[t] * n_j * b_k;
    }

  } else {
    for (int t = 0; t < nTerms; ++t) {
      int ii = terms[t].i;
      int jj = 2 + terms[t].j;
      int kk = 5 + terms[t].k;

      GDouble n_j = GPU_UVARS(jj);  // kN0..kN2
      GDouble b_k = GPU_UVARS(kk);  // kB0..kB5
      
      // if (threadIdx.x == 0 && blockIdx.x == 0) {
      //   printf(
      //       "t=%d ii=%d jj=%d kk=%d\n g_i=%g n_j=%g b_k=%g coeffs=%g\n "
      //       "kN0=%g kN1=%g kN2=%g\n "
      //       "kB0=%g kB1=%g kB2=%g kB3=%g kB4=%g kB5=%g\n",
      //       t, ii, jj, kk,
      //       g_func[ii], n_j, b_k, coeffs[t],
      //       GPU_UVARS(2), GPU_UVARS(3), GPU_UVARS(4),
      //       GPU_UVARS(5), GPU_UVARS(6), GPU_UVARS(7),
      //       GPU_UVARS(8), GPU_UVARS(9), GPU_UVARS(10)
      //   );
      // }
      H[ii] += coeffs[t] * g_func[ii] * n_j * b_k;
    }
  }

  // ----- optional Schilling–Wolf SDME -----
  GDouble W[3][2] = {
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0}
  };

  if (schillingIncluded) {
    GDouble rho000  = rhos[0];
    GDouble rho100  = rhos[1];
    GDouble rho1m10 = rhos[2];
    GDouble rho111  = rhos[3];
    GDouble rho001  = rhos[4];
    GDouble rho101  = rhos[5];
    GDouble rho1m11 = rhos[6];
    GDouble rho102  = rhos[7];
    GDouble rho1m12 = rhos[8];

    W[0][0] += 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*GPU_UVARS(6) - G_SQRT(2.0)*rho100*GPU_UVARS(7) - rho1m10*GPU_UVARS(9);
    W[1][0] += (rho111*(1. - GPU_UVARS(6)) + rho001*GPU_UVARS(6) - G_SQRT(2.0)*rho101*GPU_UVARS(7) - rho1m11*GPU_UVARS(9));
    W[2][1] += (G_SQRT(2.0)*rho102*GPU_UVARS(8) + rho1m12*GPU_UVARS(10));

    if (model == 2) {
      GDouble rho100_neg  = rhos[9];
      GDouble rho1m10_neg = rhos[10];

      GDouble rho101_neg  = rhos[11];
      GDouble rho1m11_neg = rhos[12];

      GDouble rho112_neg  = rhos[13];
      GDouble rho002_neg  = rhos[14];
      GDouble rho102_neg  = rhos[15];
      GDouble rho1m12_neg = rhos[16];

      W[0][1] += G_SQRT(2.0) * rho100_neg * GPU_UVARS(8) + rho1m10_neg * GPU_UVARS(10);
      W[1][1] += G_SQRT(2.0) * rho101_neg * GPU_UVARS(8) + rho1m11_neg * GPU_UVARS(10);
      W[2][1] += rho112_neg*(1. - GPU_UVARS(6)) + rho002_neg*GPU_UVARS(6) - G_SQRT(2.0)*rho102_neg*GPU_UVARS(7) - rho1m12_neg*GPU_UVARS(9);
    }
  }
  else{
    W[0][0] = 1.0;
  }

  GDouble W_normalization = 3.0 / (4.0 * PI);
  GDouble H_normalization = 1.0 / (4.0 * PI);

  GDouble I_total = 0.0;
  if (model != 2) {
    I_total += W_normalization * (W[0][0] + g_func[1] * W[1][0] + g_func[2] * W[2][0])
             + alpha * H_normalization * (H[0] + H[1] + H[2]);
  } else if (model == 2) {
    I_total += W_normalization * H_normalization * (
                              W[0][0] * FH[0][0] + W[0][1] * FH[0][1] +
                 g_func[1] * (W[1][0] * FH[1][0] + W[1][1] * FH[1][1]) +
                 g_func[2] * (W[2][0] * FH[2][1] + W[2][1] * FH[2][0])
               );
  }

  WCUComplex amp = { G_SQRT(G_FABS(I_total)), 0 };
  pcDevAmp[iEvent] = amp;
}

void
GPUKStarHyperonExtended_exec(
  dim3 dimGrid,
  dim3 dimBlock,
  GPU_AMP_PROTO,
  int model,
  GDouble alpha,
  GDouble* coeffs,          // host array (size nCoeffs)
  term_ijk* terms,          // host array (size nTerms)
  int nTerms,
  bool schillingIncluded,
  GDouble* rhos,            // host array size 9 if schillingIncluded else nullptr
  GDouble polAngle
) {
  // allocate device buffers
  GDouble* d_coeffs = nullptr;
  term_ijk* d_terms = nullptr;
  GDouble* d_rhos   = nullptr;

  // coeff count is NOT the same as nTerms for model==2
  const int nCoeff = (model >= 3 ? 39 : nTerms);        
  // allocate/copy coeffs (always)
  cudaMalloc((void**)&d_coeffs, nCoeff * sizeof(GDouble));              
  cudaMemcpy(d_coeffs, coeffs,  nCoeff * sizeof(GDouble), cudaMemcpyHostToDevice); 

   // allocate/copy terms only if needed
   if (nTerms > 0) {                                    
    cudaMalloc((void**)&d_terms,  nTerms * sizeof(term_ijk));
    cudaMemcpy(d_terms,  terms,   nTerms * sizeof(term_ijk), cudaMemcpyHostToDevice);
  } else {
    d_terms = nullptr;                                 
  }

  if (schillingIncluded) {
    cudaMalloc((void**)&d_rhos, 9 * sizeof(GDouble));
    cudaMemcpy(d_rhos, &rhos[0], 9 * sizeof(GDouble), cudaMemcpyHostToDevice);
  }

  // launch kernel with device pointers
  GPUKStarHyperonExtended_kernel<<< dimGrid, dimBlock >>>(
    GPU_AMP_ARGS, model, alpha, d_coeffs, d_terms, nTerms, schillingIncluded, d_rhos, polAngle
  );

  cudaDeviceSynchronize();
  // auto e = cudaGetLastError();
  // if (e != cudaSuccess) printf("CUDA kernel launch error: %s\n", cudaGetErrorString(e));

  // free device memory
  cudaFree(d_coeffs);
  if (d_terms) cudaFree(d_terms);                      
  if (d_rhos)  cudaFree(d_rhos);
}
