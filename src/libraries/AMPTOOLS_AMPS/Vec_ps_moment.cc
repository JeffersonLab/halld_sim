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
  // count the number of non-moment arguments. Sometimes moments are passed within
  // brackets "[]", so this accounts for that case.
  m_nonMomentArgs = 0; 
  for (const auto& arg : args) {
    if (arg.size() < 2 || (arg[0] != 'H' && arg[1] != 'H')) {
      m_nonMomentArgs++;
    }
  }

  m_numberOfMoments = args.size() - m_nonMomentArgs; // calculate the number of moments based on non-moment arguments
  m_moments.reserve(m_numberOfMoments);

  // Three possibilities to initialize the set of polarized moment parameters, with the following labels:
  // <alpha>: indexes the unpolarized [0], polarized real [1], and polarized imaginary moments [2],
  // <Jv>: moment quantum number originating from decay of vector meson,
  // <Lambda>: moment quantum number originating from omega helicity
  // <J>: moment quantum number originating from spin of resonance decay
  // <M>: moment quantum number originating from spin projection of resonance decay
  
  // NOTE that in each case below, the non-moment parameters are given in a specific
  // order, and always listed before the moments.
  
  // 1: Only moments are given, no polarization information
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <Moments...>
  if(m_nonMomentArgs == 0) {
    m_polInTree = true;
  }
  // 2: The beam polarization angle and fraction are fixed by the user (2 arguments)
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <polAngle> <polFraction> <Moments...> 
  else if (m_nonMomentArgs == 2) {
    m_polInTree = false;
    m_polAngle = stof( args[0].c_str() );
    m_polFraction = stof( args[1].c_str() );
  }
  // 3: Polarization read from histogram <hist> in file <rootFile> (3 arguments)
  //    Usage: amplitude <reaction>::<sum>::<ampName> Vec_ps_moment <polAngle> <rootFile> <hist> <Moments...> 
  else if(m_nonMomentArgs == 3) {
    m_polInTree = false;
    m_polAngle = stof( args[0].c_str() );
    m_polFraction = 0.; 
    TFile* f = new TFile( args[1].c_str() );
    m_polFrac_vs_E = (TH1D*)f->Get( args[2].c_str() );
    assert( m_polFrac_vs_E != NULL );
    delete f; 
  } 
  else {
    cout << " ERROR: Invalid number of arguments for Vec_ps_moment constructor" << endl;
    assert(0);
  }

  for(int i = m_nonMomentArgs; i < m_numberOfMoments + m_nonMomentArgs; i++) { 
    moment mom;

    mom.name = args[i];
    mom.H = AmpParameter( mom.name );     
    registerParameter( mom.H ); // register moment as a free parameter

    // remove brackets from moment name if they exist to parse it
    if (mom.name.front() == '[' && mom.name.back() == ']') {
      mom.name = mom.name.substr(1, mom.name.size() - 2);
    }

    // Check if the moment name has the expected length of 7 characters
    if (mom.name.length() != 7) {
      cout << " ERROR: Invalid moment name length for " << mom.name << endl;
      assert(0);
    }

    // parse the moment name to get the quantum numbers. Assumes moment name is of the
    // form "H<alpha>_<Jv><Lambda><J><M>"
    mom.alpha = stoi(mom.name.substr(1,1).data());
    mom.Jv = stoi(mom.name.substr(3,1).data());
    mom.Lambda = stoi(mom.name.substr(4,1).data());
    mom.J = stoi(mom.name.substr(5,1).data());
    mom.M = stoi(mom.name.substr(6,1).data());

    m_moments.push_back( mom ); // add the moment to the list
  }   

  // default is 3-body vector decay (omega->3pi)
  m_3pi = true; 
}

void
Vec_ps_moment::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector beam;
  TVector3 eps;
  GDouble beam_polFraction;

  if(m_polInTree){
    beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // beam polarization vector;

    beam_polFraction = eps.Mag();    
  }
  else {
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    eps.SetXYZ(cos(m_polAngle*TMath::DegToRad()), sin(m_polAngle*TMath::DegToRad()), 0.0); // beam polarization vector

    if(m_polFraction > 0.) { // for fitting with fixed polarization
      beam_polFraction = m_polFraction;
    } 
    else { // for fitting with polarization vs E_gamma from input histogram 
      int bin = m_polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
      if (bin == 0 || bin > m_polFrac_vs_E->GetXaxis()->GetNbins()) {
        beam_polFraction = 0.;
      } 
      else {
        beam_polFraction = m_polFrac_vs_E->GetBinContent(bin);
      }
    }
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
  
  userVars[kProdAngle] = locthetaphi[2];

  userVars[kBeamPolFraction] = beam_polFraction;
}


////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
Vec_ps_moment::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  // obtain all the angles from the userVars block
  GDouble beam_polFraction = userVars[kBeamPolFraction];
  GDouble cosTheta = userVars[kCosTheta];
  GDouble phi = userVars[kPhi];
  GDouble cosThetaH = userVars[kCosThetaH];
  GDouble phiH = userVars[kPhiH];
  GDouble prodAngle = userVars[kProdAngle];

  GDouble cos2prodAngle = cos(2*prodAngle);
  GDouble sin2prodAngle = sin(2*prodAngle);
  GDouble theta = acos(cosTheta) * 180./TMath::Pi(); 
  GDouble thetaH = acos(cosThetaH) * 180./TMath::Pi(); 

  GDouble total = 0; // initialize the total "amplitude" to zero
  for(int i = 0; i < m_numberOfMoments; i++) {
    moment mom = m_moments[i];
    // extract the quantum numbers
    int alpha = mom.alpha;
    int Jv = mom.Jv;
    int Lambda = mom.Lambda;
    int J = mom.J;
    int M = mom.M;

    // initialize angular info with the common factor
    GDouble angle = (2*J + 1) / (4*TMath::Pi()) * (2*Jv + 1) / (4*TMath::Pi());

    // calculate the angular distributions depending on the polarization component
    // NOTE: the multiplication of two Wigner D functions gives 2 wignerDsmall functions,
    // and an exponential. We enforce that the H0, H1 (H2) moments are purely real 
    // (imaginary) by taking the real (imaginary) part of that exponential.
    if(alpha == 0) { // unpolarized component
      angle *= (
        wignerDSmall( J, M, Lambda, theta ) * 
        wignerDSmall( Jv, Lambda, 0, thetaH ) * 
        cos( M*phi + Lambda*phiH )
      );
    }
    else if(alpha == 1) {
      angle *= -1.0 * beam_polFraction * cos2prodAngle * 
        wignerDSmall( J, M, Lambda, theta ) * 
        wignerDSmall( Jv, Lambda, 0, thetaH ) * 
        cos( M*phi + Lambda*phiH );
    }
    else if(alpha == 2) {
      angle *= -1.0 * beam_polFraction * sin2prodAngle * 
        wignerDSmall( J, M, Lambda, theta ) * 
        wignerDSmall( Jv, Lambda, 0, thetaH ) * 
        sin( M*phi + Lambda*phiH );
    }

    total += angle * mom.H;
  }

  // since AmpTools is hard coded to handle squared amplitudes, we return the square root of the total
  return complex< GDouble >( sqrt(fabs(total)) );
}

#ifdef GPU_ACCELERATION

void
Vec_ps_moment::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  // convert vector to array for GPU
  GDouble H[m_numberOfMoments];
  moment moments[m_numberOfMoments];
  for(int i=0; i<m_numberOfMoments; i++)
  {
    H[i] = m_moments[i].H;
    moments[i] = m_moments[i];
  }  

  GPUVec_ps_moment_exec( dimGrid, dimBlock, GPU_AMP_ARGS, H, moments, m_numberOfMoments );

}

#endif


