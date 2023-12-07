#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPSMoment.h"
#include "AMPTOOLS_AMPS/wignerD.h"

#include "TFile.h"

TwoPSMoment::TwoPSMoment( const vector< string >& args ) :
   UserAmplitude< TwoPSMoment >( args )
{
   m_maxL = atoi( args[0].c_str() );
   m_nMoments = atoi( args[1].c_str() );
   for(int imom = 0; imom < m_nMoments; imom++) { 
	   m_H.push_back( AmpParameter( args[imom+2] ) );
	   
	   string name = m_H[imom].name();
	   
	   int index = atoi(name.substr(1,1).data()) * 100;
	   index += atoi(name.substr(3,1).data()) * 10;
	   index += atoi(name.substr(4,1).data()); 
	   	   
	   // global moment indices
	   m_indices.push_back( index );
   }

   for(int imom = 0; imom < m_nMoments; imom++) 
	   registerParameter( m_H[imom] );

   // Three possibilities to initialize this amplitude:
   // (with <alpha>, <l>, <m> defining the polarized moments)
   //
   // 1: three arguments, polarization information must be included in beam photon four vector
   //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPSMoment <maxL> <nMoments> <Moments...>
   if(args.size() == size_t(m_nMoments + 2)) {
      m_polInTree = true;
   
   // 2: five arguments, polarization fixed per amplitude and passed as flag
   //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPSMoment <maxL> <nMoments> <Moments...> <polAngle> <polFraction>
   } else if(args.size() == size_t(m_nMoments + 4)) {
      m_polInTree = false;
      m_polAngle = atof( args[m_nMoments+2].c_str() );
      m_polFraction = atof( args[m_nMoments+3].c_str() );
   
   // 3: eight arguments, read polarization from histogram <hist> in file <rootFile>
   //    Usage: amplitude <reaction>::<sum>::<ampName> TwoPSMoment <maxL> <nMoments> <Moments...> <polAngle> <polFraction=0.> <rootFile> <hist>
   } else if(args.size() == size_t(m_nMoments + 6)) {
      m_polInTree = false;
      m_polAngle = atof( args[m_nMoments+2].c_str() );
      m_polFraction = 0.; 
      TFile* f = new TFile( args[m_nMoments+4].c_str() );
      m_polFrac_vs_E = (TH1D*)f->Get( args[m_nMoments+5].c_str() );
      assert( m_polFrac_vs_E != NULL );
   }
}


complex< GDouble >
TwoPSMoment::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

   GDouble pGamma = userVars[kPgamma];
   GDouble cosTheta = userVars[kCosTheta];
   GDouble phi = userVars[kPhi];
   GDouble bigPhi = userVars[kBigPhi];

   GDouble cos2bigPhi = cos(2*bigPhi);
   GDouble sin2bigPhi = sin(2*bigPhi);
   GDouble theta = acos(cosTheta) * 180./TMath::Pi(); 

   // compute required wignerDSmall values
   GDouble wigner[10*(m_maxL+1)];
   for(int iL = 0; iL <= m_maxL; iL++) {
	   for(int iM = 0; iM <= iL; iM++) {
		   wigner[iL*10 + iM] = wignerDSmall( iL, iM, 0, theta );
	   }
   }

   GDouble total = 0;
   for(int imom = 0; imom < m_nMoments; imom++) { 
	   	   
	   int alpha = m_indices[imom] / 100;
	   int L = m_indices[imom] / 10 % 10;
	   int M = m_indices[imom] % 10;
	   int LM = m_indices[imom] % 100;

	   GDouble mom = 0;
#if 0
	   mom = 2.0 * sqrt( (2*L + 1) / (4*TMath::Pi()) );
	   if(alpha == 0)
		   mom *= sphericalHarmonics[L][M].real();
	   else if(alpha == 1) 
		   mom *= pGamma * cos2bigPhi * sphericalHarmonics[L][M].real();
	   else if(alpha == 2) 
		   mom *= -1.0 * pGamma * sin2bigPhi * sphericalHarmonics[L][M].imag();	   
#endif

#if 1
	   mom = 2.0 * (2*L + 1) / (4*TMath::Pi());
	   if(alpha == 0)
		   mom *= wigner[LM] * cos( M*phi );
	   else if(alpha == 1) 
		   mom *= pGamma * cos2bigPhi * wigner[LM] * cos( M*phi );
	   else if(alpha == 2) 
		   mom *= -1.0 * pGamma * sin2bigPhi * wigner[LM] * sin( M*phi );
#endif
	   
	   // m = 0 only non-zero for alpha = 0, 1 but half the size of other m-projections
	   if(M == 0 && alpha < 2) mom *= 0.5;
	   
	   total += m_H[imom]*mom;
   }   
	   
   return complex< GDouble >( sqrt(fabs(total)) );
}

void
TwoPSMoment::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
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
TwoPSMoment::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

   // convert vector to array for GPU
   GDouble H[m_nMoments];
   int indices[m_nMoments];
   for(int i=0; i<m_nMoments; i++){
      H[i] = m_H[i];
      indices[i] = m_indices[i];
   }

   GPUTwoPSMoment_exec( dimGrid, dimBlock, GPU_AMP_ARGS, H, indices, m_nMoments, m_maxL );
}
#endif
