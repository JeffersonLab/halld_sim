#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/VecPSMoment.h"
#include "AMPTOOLS_AMPS/wignerD.h"

#include "TFile.h"

VecPSMoment::VecPSMoment( const vector< string >& args ) :
   UserAmplitude< VecPSMoment >( args )
{
   m_maxJ = atoi( args[0].c_str() );
   m_nMoments = atoi( args[1].c_str() );
   for(int imom = 0; imom < m_nMoments; imom++) { 
	   m_H.push_back( AmpParameter( args[imom+2] ) );
	   
	   string name = m_H[imom].name();
	   m_alpha.push_back( atoi(name.substr(1,1).data()) );
	   m_S.push_back( atoi(name.substr(3,1).data()) );
	   m_Lambda.push_back( atoi(name.substr(4,1).data()) );
	   m_J.push_back( atoi(name.substr(5,1).data()) );
	   m_M.push_back( atoi(name.substr(6,1).data()) );
   }

   for(int imom = 0; imom < m_nMoments; imom++) 
	   registerParameter( m_H[imom] );

   // Three possibilities to initialize this amplitude:
   // (with <alpha>, <l>, <m> defining the polarized moments)
   //
   // 1: three arguments, polarization information must be included in beam photon four vector
   //    Usage: amplitude <reaction>::<sum>::<ampName> VecPSMoment <maxL> <nMoments> <Moments...>
   if(args.size() == size_t(m_nMoments + 2)) {
      m_polInTree = true;
   
   // 2: five arguments, polarization fixed per amplitude and passed as flag
   //    Usage: amplitude <reaction>::<sum>::<ampName> VecPSMoment <maxJ> <nMoments> <Moments...> <polAngle> <polFraction>
   } else if(args.size() == size_t(m_nMoments + 4)) {
      m_polInTree = false;
      m_polAngle = atof( args[m_nMoments+2].c_str() );
      m_polFraction = atof( args[m_nMoments+3].c_str() );
   
   // 3: eight arguments, read polarization from histogram <hist> in file <rootFile>
   //    Usage: amplitude <reaction>::<sum>::<ampName> VecPSMoment <maxJ> <nMoments> <Moments...> <polAngle> <polFraction=0.> <rootFile> <hist>
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
VecPSMoment::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

   GDouble pGamma = userVars[kPgamma];
   GDouble cosTheta = userVars[kCosTheta];
   GDouble phi = userVars[kPhi];
   GDouble cosThetaH = userVars[kCosThetaH];
   GDouble phiH = userVars[kPhiH];
   GDouble bigPhi = userVars[kBigPhi];

   GDouble theta = acos(cosTheta) * 180./TMath::Pi(); 
   GDouble thetaH = acos(cosThetaH) * 180./TMath::Pi(); 

   GDouble total = 0;
   for(int imom = 0; imom < m_nMoments; imom++) { 
	   	   
	   int alpha = m_alpha[imom];
	   int S = m_S[imom];
	   int Lambda = m_Lambda[imom];
	   int J = m_J[imom];
	   int M = m_M[imom];

	   GDouble mom = sqrt( (2*J + 1) / (4*TMath::Pi() ) ) * sqrt( (2*S + 1) / (4*TMath::Pi() ) );
	   if(alpha == 0)
		   mom *= wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( S, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );
	   else if(alpha == 1) 
		   mom *= -1.0 * pGamma * cos(2*bigPhi) * wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( S, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );
	   else if(alpha == 2) 
		   mom *= -1.0 * pGamma * sin(2*bigPhi) * wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( S, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );
	   	   
	   total += mom*m_H[imom];
   }   

   return complex< GDouble >( sqrt(fabs(total)) );
}

void
VecPSMoment::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
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

   // place holders
   userVars[kCosThetaH] = 0;
   userVars[kPhiH] = 0;

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
VecPSMoment::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

   // convert vector to array for GPU
   GDouble H[m_nMoments];
   int alpha[m_nMoments];
   int S[m_nMoments];
   int Lambda[m_nMoments];
   int J[m_nMoments];
   int M[m_nMoments];
   for(int i=0; i<m_nMoments; i++){
      H[i] = m_H[i];
      alpha[i] = m_alpha[i];
      S[i] = m_S[i];
      Lambda[i] = m_Lambda[i];
      J[i] = m_J[i];
      M[i] = m_M[i];
   }


   GPUVecPSMoment_exec( dimGrid, dimBlock, GPU_AMP_ARGS, H, alpha, S, Lambda, J, M, m_nMoments );
}
#endif
