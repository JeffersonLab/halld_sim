
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/VecPs_3pi_refl.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/vecPsAngles.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"

#include "UTILITIES/BeamProperties.h"

// Function to check if a string is a valid number
static bool isValidNumber(const string& argInput, double &value){
    char* end = nullptr;
    errno = 0;  // reset global error
    value = std::strtod(argInput.c_str(), &end);

    // Check if 
    // (1) no conversion was performed 
    // (2) there are leftover characters
    // (3) an overflow/underflow occurred   
    if(end == argInput.c_str() || *end != '\0' || errno != 0) {
        return false;  // not a valid number
    }
    // If end points to the end of string, it's fully numeric
    return true;
}

static double parseValidatedNumber(const string& label, const string& argInput){
    double tmpValue = 0.0;
    if(!isValidNumber(argInput, tmpValue)){
      throw std::invalid_argument("Vec_ps_refl: invalid " + label + ": " + argInput);
    }
    return tmpValue;
}


VecPs_3pi_refl::VecPs_3pi_refl( const vector< string >& args ) :
UserAmplitude< VecPs_3pi_refl >( args ){

  // 5 possibilities to initialize this amplitude:
  // <J>: total spin, 
  // <m>: spin projection, 
  // <l>: partial wave, 
  // <r>: +1/-1 for real/imaginary part; 
  // <s>: +1/-1 sign in P_gamma term
  
  // Resonance spin J
  m_j = static_cast<int>(parseValidatedNumber("J", args[0])); 
  // Spin projection (Lambda)
  m_m = static_cast<int>(parseValidatedNumber("M", args[1])); 
  // Partial wave L
  m_l = static_cast<int>(parseValidatedNumber("L", args[2])); 
  // Real (+1) or imaginary (-1)
  m_r = static_cast<int>(parseValidatedNumber("Re/Im", args[3])); 
  // Sign for polarization in amplitude
  m_s = static_cast<int>(parseValidatedNumber("P_gamma sign", args[4])); 

  // make sure values are reasonable
  assert( abs( m_m ) <= m_j );
  // m_r = +1 for real
  // m_r = -1 for imag
  assert( abs( m_r ) == 1 );
  // m_s = +1 for 1 + Pgamma
  // m_s = -1 for 1 - Pgamma
  assert( abs( m_s ) == 1 );

  // Default polarization information stored in tree
  m_polInTree = true;

  // Default is 2-body vector decay (set flag in config file for omega->3pi)
  m_3pi = false;

  // Loop over any additional amplitude arguments to change defaults
  for(uint ioption=5; ioption<args.size(); ioption++) {
	  TString option = args[ioption].c_str();
    // Polarization provided in configuration file
    if(ioption==5){
      m_polInTree = false;

      polAngle = parseValidatedNumber("polarization angle", args[5]);    

      TString polOption = args[6].c_str();
      if(polOption.Contains(".root")){
        polFraction = 0.0;
        TFile* f = new TFile(polOption);
        polFrac_vs_E = (TH1D*)f->Get(args[7].c_str());
        if(polFrac_vs_E  != nullptr ){
          throw std::runtime_error(
            "VecPs_3pi_refl ERROR: Could not find histogram '" + args[7] +
            "' in file " + std::string(polOption.Data()));
        }
      }
      else{
        polFraction = parseValidatedNumber("polarization fraction", args[6]);
      }
    }
    // Check for omega->3pi option
	  if(option.EqualTo("omega3pi")) m_3pi = true;
  }  
}

void
VecPs_3pi_refl::calcUserVars( GDouble** pKin, GDouble* userVars ) const{

  TLorentzVector beam;
  TVector3 eps;
  GDouble beam_polFraction;
  GDouble beam_polAngle;

  if(m_polInTree){
    beam.SetPxPyPzE( 0.0, 0.0, pKin[0][0], pKin[0][0]);
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.0); // beam polarization vector;

    beam_polFraction = eps.Mag();
    beam_polAngle = eps.Phi();
  }
  else{
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    beam_polAngle = polAngle;

    if(polFraction > 0.0){ // for fitting with fixed polarization
	    beam_polFraction = polFraction;
    }
    else{ // for fitting with polarization vs E_gamma from input histogram 
	    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
	    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		    beam_polFraction = 0.0;
	    } else
		    beam_polFraction = polFrac_vs_E->GetBinContent(bin);
    }
  }
  
  TLorentzVector recoil( pKin[1][1]+pkin[5][1], pKin[1][2]+pkin[5][2], pKin[1][3]+pkin[5][3], pKin[1][0]+pKin[5][0] ); 

  // Fill in four-vectors for final state particles
  // 1st after proton is always the pseudoscalar meson
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); 
  // Compute vector meson from its decay products
  // Make sure the order of daughters is correct in the config file!
  TLorentzVector vec, vec_daught1, vec_daught2; 


  if(m_3pi){
    // Omega ps proton, omega -> 3pi (6 particles)
    // Omega pi- Delta++, omega -> 3pi (7 particles)
	  TLorentzVector pi0(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  TLorentzVector pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  TLorentzVector pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
	  vec = pi0 + pip + pim;
	  vec_daught1 = pip;
	  vec_daught2 = pim;
  }
  else{
	  // Omega ps proton, omega -> pi0 g (4 particles)
	  // Omega pi- Delta++, omega -> pi0 g (5 particles)

	  // (vec 2-body) ps proton, vec 2-body -> pipi, KK (5 particles)
	  // (vec 2-body) pi- Delta++, vec 2-body -> pipi, KK (6 particles)
	  // (vec 2-body) K+ Lambda, vec 2-body -> Kpi (6 particles)
	  vec_daught1 = TLorentzVector(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  vec_daught2 = TLorentzVector(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  vec = vec_daught1 + vec_daught2;
  }

  // Final meson system P4
  TLorentzVector X = vec + ps;

  ///////////////// Boost Particles and Get Angles/////////////////////

  TLorentzVector target(0,0,0,0.938);
  // Helicity coordinate system
  TLorentzVector beamP = beam + target;

  // Calculate decay angles for X in helicity frame (same for all vectors)
  // Change getXDecayAngles to get Gottfried-Jackson angles if needed
  // Note: it also calculates the production angle
  vector <double> xDecayAngles = getXDecayAngles( beam_polAngle, beam, beamP, X, vec);

  // Calculate vector decay angles (unique for each vector)
  vector <double> vectorDecayAngles;
  if(m_3pi){
    vectorDecayAngles = getVectorDecayAngles( beamP, X, vec,
                                              vec_daught1, vec_daught2);
  }
  else{
    vectorDecayAngles = getVectorDecayAngles( beamP, X, vec,
                                        vec_daught1, TLorentzVector(0,0,0,0));
  }

  GDouble cosTheta = TMath::Cos(xDecayAngles[0]);
  GDouble phi = xDecayAngles[1];
  GDouble prod_angle = xDecayAngles[2]; // bigPhi
  GDouble cosThetaH = TMath::Cos(vectorDecayAngles[0]);
  GDouble phiH = vectorDecayAngles[1];
  GDouble m_X = X.M();
  GDouble m_vec = vec.M();
  GDouble m_ps = ps.M();

  complex <GDouble> amplitude(0,0);
  complex <GDouble> i(0,1);

  for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
	  GDouble hel_amp = clebschGordan(m_l, 1, 0, lambda, m_j, lambda);
          amplitude += conj(wignerD(m_j, m_m, lambda, cosTheta, phi)) *
                       hel_amp * conj(wignerD(1, lambda, 0, cosThetaH, phiH));
  }

  GDouble factor = sqrt(1 + m_s * beam_polFraction);
  complex <GDouble> zjm = 0;
  // - -> + in prod_angle
  complex <GDouble> rotateY = polar((GDouble)1., (GDouble)(-1. * prod_angle ));  

  if (m_r == 1)
	  zjm = real(amplitude * rotateY);
  if (m_r == -1) 
	  zjm = i*imag(amplitude * rotateY);

  // E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
  double kinFactor = barrierFactor(m_X, m_l, m_vec, m_ps);
  //kinFactor *= sqrt(3.) * sqrt(2.*m_l + 1.);
  factor *= kinFactor;

  userVars[uv_ampRe] = ( factor * zjm ).real();
  userVars[uv_ampIm] = ( factor * zjm ).imag();

  return;
}


/////////////////////// Amplitude Calculation //////////////////////////

complex< GDouble >
VecPs_3pi_refl::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  return complex< GDouble >( userVars[uv_ampRe], userVars[uv_ampIm] );
}

void VecPs_3pi_refl::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
VecPs_3pi_refl::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUVecPs_3pi_refl_exec( dimGrid, dimBlock, GPU_AMP_ARGS );

}

#endif


