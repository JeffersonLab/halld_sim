#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Zlm.h"
#include "AMPTOOLS_AMPS/wignerD.h"

#include "TFile.h"

Zlm::Zlm( const vector< string >& args ) :
   UserAmplitude< Zlm >( args )
{
   assert( args.size() == 4 || args.size() == 6 || args.size() == 8 );

   m_j = atoi( args[0].c_str() );
   m_m = atoi( args[1].c_str() );
   m_r = atoi( args[2].c_str() );
   m_s = atoi( args[3].c_str() );

   // Three possibilities to initialize this amplitude:
   // (with <J>: total spin, <m>: spin projection, <r>: +1/-1 for real/imaginary part; <s>: +1/-1 sign in P_gamma term)
   //
   // 1: four arguments, polarization information must be included in beam photon four vector
   //    Usage: amplitude <reaction>::<sum>::<ampName> Zlm <J> <m> <r> <s>
   if(args.size() == 4) {
      m_polInTree = true;
   
   // 2: six arguments, polarization fixed per amplitude and passed as flag
   //    Usage: amplitude <reaction>::<sum>::<ampName> Zlm <J> <m> <r> <s> <polAngle> <polFraction>
   } else if(args.size() == 6) {
      m_polInTree = false;
      m_polAngle = atof( args[4].c_str() );
      m_polFraction = atof( args[5].c_str() );
   
   // 3: eight arguments, read polarization from histogram <hist> in file <rootFile>
   //    Usage: amplitude <reaction>::<sum>::<ampName> Zlm <J> <m> <r> <s> <polAngle> <polFraction=0.> <rootFile> <hist>
   } else if(args.size() == 8) {
      m_polInTree = false;
      m_polAngle = atof( args[4].c_str() );
      m_polFraction = 0.; 
      TFile* f = new TFile( args[6].c_str() );
      m_polFrac_vs_E = (TH1D*)f->Get( args[7].c_str() );
      assert( m_polFrac_vs_E != NULL );
   }

   // make sure values are reasonable
   assert( abs( m_m ) <= m_j );
   // m_r = +1 for real
   // m_r = -1 for imag
   assert( abs( m_r ) == 1 );
   // m_s = +1 for 1 + Pgamma
   // m_s = -1 for 1 - Pgamma
   assert( abs( m_s ) == 1 );
}


complex< GDouble >
Zlm::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

   GDouble pGamma = userVars[kPgamma];
   GDouble cosTheta = userVars[kCosTheta];
   GDouble phi = userVars[kPhi];
   GDouble bigPhi = userVars[kBigPhi];

   GDouble factor = sqrt(1 + m_s * pGamma);
   GDouble zlm = 0;
   complex< GDouble > rotateY = polar(1., -1.*bigPhi);
   if (m_r == 1)
      zlm = real(Y( m_j, m_m, cosTheta, phi ) * rotateY);
   if (m_r == -1)
      zlm = imag(Y( m_j, m_m, cosTheta, phi ) * rotateY);

   return complex< GDouble >( factor * zlm );
}

void
Zlm::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
   TLorentzVector beam;
   TVector3 eps;

   if(m_polInTree) {
      beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
      eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // makes default output gen_amp trees readable as well (without transforming)
   } else {
      beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
      eps.SetXYZ(cos(m_polAngle*TMath::DegToRad()), sin(m_polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
   }

   TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
   TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
   TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 

   TLorentzVector resonance = p1 + p2;

   TLorentzRotation resRestBoost( -resonance.BoostVector() );

   TLorentzVector beam_res   = resRestBoost * beam;
   TLorentzVector recoil_res = resRestBoost * recoil;
   TLorentzVector p1_res = resRestBoost * p1;

   // Helicity frame
   TVector3 z = -1. * recoil_res.Vect().Unit();
   // or GJ frame?
   // TVector3 z = beam_res.Vect().Unit();

   // normal to the production plane
   TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

   TVector3 x = y.Cross(z);

   TVector3 angles( (p1_res.Vect()).Dot(x),
         (p1_res.Vect()).Dot(y),
         (p1_res.Vect()).Dot(z) );

   userVars[kCosTheta] = angles.CosTheta();
   userVars[kPhi] = angles.Phi();

   userVars[kBigPhi] = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));


   GDouble pGamma;
   if(m_polInTree) {
      pGamma = eps.Mag();
   } else {
      if(m_polFraction > 0.) { // for fitting with constant polarization 
         pGamma = m_polFraction;
      } else{
         int bin = m_polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
         if (bin == 0 || bin > m_polFrac_vs_E->GetXaxis()->GetNbins()){
            pGamma = 0.;
         } else 
            pGamma = m_polFrac_vs_E->GetBinContent(bin);
      }
   }

   userVars[kPgamma] = pGamma;
}

#ifdef GPU_ACCELERATION
void
Zlm::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

   GPUZlm_exec( dimGrid, dimBlock, GPU_AMP_ARGS, m_j, m_m, m_r, m_s );
}
#endif
