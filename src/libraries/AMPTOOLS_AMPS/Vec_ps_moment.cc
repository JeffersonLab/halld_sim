/* Simplistic model for fitting vector-pseudoscalar data with linearly polarized moments 

The file currently assumes the vector is an omega decaying to 3pi, and so uses the
m_3pi flag, just as is done in vec_ps_refl.cc.

Original author: Edmundo Barriga
Modified by: Kevin Scheuer
*/
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Vec_ps_moment.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/omegapiAngles.h"

Vec_ps_moment::Vec_ps_moment( const vector< string >& args ) :
  UserAmplitude< Vec_ps_moment >( args )
{
  m_maxJ = atoi( args[0].c_str() ); // explicitly set the maximum J value TODO: why?
  m_nMoments = atoi( args[1].c_str() );  // number of moment parameters

  for(int imom = 0; imom < m_nMoments; imom++) { 
    string parameter = args[imom+2];
    string name = parameter.name();
    m_H.push_back( AmpParameter( parameter ) );

    // TODO: why is this indexing needed?
    int index = atoi(name.substr(1,1).data()) * 10000;
    index += atoi(name.substr(3,1).data()) * 1000;
    index += atoi(name.substr(4,1).data()) * 100;
    index += atoi(name.substr(5,1).data()) * 10;
    index += atoi(name.substr(6,1).data());

    // global moment indices
    m_indices.push_back( index );

    // register free parameters so the framework knows about them
    registerParameter( parameter );
  }   

  // Three possibilities to initialize the set of polarized moment parameters, with the following labels:
  // <alpha>: indexes the unpolarized [0], polarized real [1], and polarized imaginary moments [2],
  // <Jv>: moment quantum number originating from decay of vector meson,
  // <Lambda>: moment quantum number originating from omega helicity
  // <J>: moment quantum number originating from spin of resonance decay
  // <M>: moment quantum number originating from spin projection of resonance decay

  // note that below the vector of moments is listed as a single argument, but in reality
  // it is a list of arguments, one for each moment parameter i.e. H<alpha>_<Jv><Lambda><J><M>
  
  // 1: three arguments, polarization information must be included in a beam photon four vector
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <maxJ> <nMoments> <Moments...>
  if(args.size() == size_t(m_nMoments + 2)) {
    m_polInTree = true;
  }
  // 2: five arguments, fixed polarization defined by user
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <maxJ> <nMoments> <Moments...> <polAngle> <polFraction>
  else if(args.size() == size_t(m_nMoments + 4)) {
    m_polInTree = false;
    m_polAngle = atof( args[m_nMoments+2].c_str() );
    m_polFraction = atof( args[m_nMoments+3].c_str() );  
  }
  // 3: seven arguments, read polarization from histogram <hist> in file <rootFile>
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <maxJ> <nMoments> <Moments...> <polAngle> <polFraction.> <rootFile> <hist>  
  else if(args.size() == size_t(m_nMoments + 6)) {
    m_polInTree = false;
    m_polAngle = atof( args[m_nMoments+2].c_str() );
    m_polFraction = 0.; 
    TFile* f = new TFile( args[m_nMoments+4].c_str() );
    m_polFrac_vs_E = (TH1D*)f->Get( args[m_nMoments+5].c_str() );
    assert( m_polFrac_vs_E != NULL );
  }

  // default is 3-body vector decay (omega->3pi)
  m_3pi = true; 
}

void
Vec_ps_moment::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector beam;
  TVector3 eps;

  if(m_polInTree){
    beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // beam polarization vector;
  }
  else {
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    eps.SetXYZ(cos(m_polAngle*TMath::DegToRad()), sin(m_polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  }

  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 

  // common vector and pseudoscalar P4s
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); // 1st after proton
  TLorentzVector vec, vec_daught1, vec_daught2; // compute for each final state below 

  // omega ps proton, omega -> 3pi (6 particles)
  // omega pi- Delta++, omega -> 3pi (7 particles)
  if(m_3pi) {
    TLorentzVector pi0(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
    TLorentzVector pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
    TLorentzVector pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
    vec = pi0 + pip + pim;
    vec_daught1 = pip;
    vec_daught2 = pim;
  }
  else {
    // omega ps proton, omega -> pi0 g (4 particles)
    // omega pi- Delta++, omega -> pi0 g (5 particles)

    // (vec 2-body) ps proton, vec 2-body -> pipi, KK (5 particles)
    // (vec 2-body) pi- Delta++, vec 2-body -> pipi, KK (6 particles)
    // (vec 2-body) K+ Lambda, vec 2-body -> Kpi (6 particles)
    vec_daught1 = TLorentzVector(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
    vec_daught2 = TLorentzVector(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
    vec = vec_daught1 + vec_daught2;
  }

  // final meson system P4
  TLorentzVector X = vec + ps;

  //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  TLorentzVector target(0,0,0,0.938);
  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;

  // Calculate decay angles in helicity frame (same for all vectors)
  vector <double> locthetaphi = getomegapiAngles(eps.Phi(), vec, X, beam, Gammap);

  // Calculate vector decay angles (unique for each vector)
  vector <double> locthetaphih;
  if(m_3pi) locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, vec_daught2);
  else locthetaphih = getomegapiAngles(vec_daught1, vec, X, Gammap, TLorentzVector(0,0,0,0));
  
  userVars[kCosTheta] = TMath::Cos(locthetaphi[0]);
  userVars[kPhi] = locthetaphi[1];
  
  userVars[kCosThetaH] = TMath::Cos(locthetaphih[0]);
  userVars[kPhiH] = locthetaphih[1];
  
  userVars[kBigPhi] = locthetaphi[2];
  
  GDouble pGamma;
  if(m_polInTree) {
    pGamma = eps.Mag();
  } 
  else {
    if(m_polFraction > 0.) { // for fitting with constant polarization 
      pGamma = m_polFraction;
    } 
    else {
      int bin = m_polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
      if (bin == 0 || bin > m_polFrac_vs_E->GetXaxis()->GetNbins()) {
        pGamma = 0.;
      } 
      else {
        pGamma = m_polFrac_vs_E->GetBinContent(bin);
      }
    }
  }

  userVars[kPgamma] = pGamma;
}

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
Vec_ps_moment::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  // obtain all the angles from the userVars block
  GDouble pGamma = userVars[kPgamma];
  GDouble cosTheta = userVars[kCosTheta];
  GDouble phi = userVars[kPhi];
  GDouble cosThetaH = userVars[kCosThetaH];
  GDouble phiH = userVars[kPhiH];
  GDouble bigPhi = userVars[kBigPhi];

  GDouble cos2bigPhi = cos(2*bigPhi);
  GDouble sin2bigPhi = sin(2*bigPhi);
  GDouble theta = acos(cosTheta) * 180./TMath::Pi(); 
  GDouble thetaH = acos(cosThetaH) * 180./TMath::Pi(); 

  GDouble total = 0; // initialize the total "amplitude" to zero
  for(int imom = 0; imom < m_nMoments; imom++)
  { 
    int alpha = m_indices[imom] / 10000;
    int Jv = m_indices[imom] / 1000 % 10;
    int Lambda = m_indices[imom] / 100 % 10;
    int J = m_indices[imom] / 10 % 10;
    int M = m_indices[imom] %10;

    // initialize moment with the common factor
    GDouble mom =(2*J + 1) / (4*TMath::Pi()) * (2*Jv + 1) / (4*TMath::Pi());

    // handle what intensity component the moment contributes to
    if(alpha == 0)
      mom *= wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( Jv, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );
    else if(alpha == 1) 
      mom *= -1.0 * pGamma * cos2bigPhi * wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( Jv, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );
    else if(alpha == 2) 
      mom *= -1.0 * pGamma * sin2bigPhi * wignerDSmall( J, M, Lambda, theta ) * wignerDSmall( Jv, Lambda, 0, thetaH ) * cos( M*phi + Lambda*phiH );

    total += mom * m_H[imom];
  } 

  return complex< GDouble >( sqrt(fabs(total)) );
}

#ifdef GPU_ACCELERATION

void
Vec_ps_moment::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  // convert vector to array for GPU
  GDouble H[m_nMoments];
  int indices[m_nMoments];
  for(int i=0; i<m_nMoments; i++)
  {
    H[i] = m_H[i];
    indices[i] = m_indices[i];
  }

  GPUVec_ps_moment_exec( dimGrid, dimBlock, GPU_AMP_ARGS, H, indices, m_nMoments );

}

#endif


