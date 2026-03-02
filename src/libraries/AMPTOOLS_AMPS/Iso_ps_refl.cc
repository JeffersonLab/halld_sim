
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Iso_ps_refl.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/decayAngles.h"
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
      throw std::invalid_argument("Iso_ps_refl: invalid " + label + ": " + argInput);
    }
    return tmpValue;
}


Iso_ps_refl::Iso_ps_refl( const vector< string >& args ) :
UserAmplitude< Iso_ps_refl >( args ){

  // This function is only for the two-body Isobar decay: xi -> pipi
  
  // 6 possibilities to initialize this amplitude:
  // <S>: isobar spin
  // <J>: total spin, 
  // <m>: total spin projection, 
  // <l>: partial wave, 
  // <r>: +1/-1 for real/imaginary part; 
  // <s>: +1/-1 sign in P_gamma term

  // Isobar spin S
  m_s = static_cast<int>(parseValidatedNumber("S",args[0]));  
  // Resonance spin J
  m_j = static_cast<int>(parseValidatedNumber("J", args[1])); 
  // Spin projection (Lambda)
  m_m = static_cast<int>(parseValidatedNumber("M", args[2])); 
  // Partial wave L
  m_l = static_cast<int>(parseValidatedNumber("L", args[3])); 
  // Real (+1) or imaginary (-1)
  m_real = static_cast<int>(parseValidatedNumber("Re/Im", args[4])); 
  // Sign for polarization in amplitude
  m_sign = static_cast<int>(parseValidatedNumber("P_gamma sign", args[5])); 

  
  // make sure values are reasonable
  assert(m_s >= 0 && m_s <= 2);
  assert( abs( m_m ) <= m_j );
  // m_real = +1 for real
  // m_real = -1 for imag
  assert( abs( m_real ) == 1 );
  // m_sign = +1 for 1 + Pgamma
  // m_sign = -1 for 1 - Pgamma
  assert( abs( m_sign ) == 1 );

  // Default polarization information stored in tree
  m_polInTree = true;

  // Loop over any additional amplitude arguments to change defaults
  for(uint ioption=6; ioption<args.size(); ioption++) {
	  TString option = args[ioption].c_str();
    // Polarization provided in configuration file
    if(ioption==6){
      m_polInTree = false;

      polAngle = parseValidatedNumber("polarization angle", args[6]);    

      TString polOption = args[7].c_str();
      if(polOption.Contains(".root")){
        polFraction = 0.0;
        TFile* f = new TFile(polOption);
        polFrac_vs_E = (TH1D*)f->Get(args[8].c_str());
        if(polFrac_vs_E  != nullptr ){
          throw std::runtime_error(
            "Iso_ps_refl ERROR: Could not find histogram '" + args[8] +
            "' in file " + std::string(polOption.Data()));
        }
      }
      else{
        polFraction = parseValidatedNumber("polarization fraction", args[7]);
      }
    }
  }  
}

void
Iso_ps_refl::calcUserVars( GDouble** pKin, GDouble* userVars ) const{

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


  // Fill in four-vectors for final state particles
  TLorentzVector target(0,0,0,0.938272);   // Fixed target
  TLorentzVector recoil( pKin[1][1]+pKin[5][1], pKin[1][2]+pKin[5][2], pKin[1][3]+pKin[5][3], pKin[1][0]+pKin[5][0] ); // Recoiling baryon
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);     // 1st after proton is always the pseudoscalar meson

  
  // Compute isobar meson from its decay products (xi -> pipi)
  // Make sure the order of daughters is correct in the config file!
  TLorentzVector iso_daught1(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector iso_daught2(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector iso = iso_daught1 + iso_daught2;

  // Final meson system P4
  TLorentzVector X = iso + ps;


  //Calculate production angle in the Gottfried-Jackson frame
  GDouble prod_angle = getPhiProd(beam_polAngle, X, beam, target, 2, true);

  // Calculate decay angles for X in the Gottfried-Jackson frame and for Isobar in the Helicity frame  
  vector <double> thetaPhiAnglesTwoStep;
  thetaPhiAnglesTwoStep = getTwoStepAngles(X, iso, iso_daught1, TLorentzVector(0,0,0,0), beam, target, 2, true);

  
  
  GDouble cosTheta = TMath::Cos(thetaPhiAnglesTwoStep[0]);
  GDouble phi = thetaPhiAnglesTwoStep[1];
  GDouble cosThetaH = TMath::Cos(thetaPhiAnglesTwoStep[2]);
  GDouble phiH = thetaPhiAnglesTwoStep[3];
  GDouble m_X = X.M();
  GDouble m_iso = iso.M();
  GDouble m_ps = ps.M();

  complex <GDouble> amplitude(0,0);
  complex <GDouble> i(0,1);

  for (int lambda = -m_s; lambda <= m_s; lambda++) { // sum over helicities
	  GDouble hel_amp = clebschGordan(m_l, m_s, 0, lambda, m_j, lambda);
          amplitude += conj(wignerD(m_j, m_m, lambda, cosTheta, phi))*hel_amp*conj(wignerD(m_s, lambda, 0, cosThetaH, phiH));
  }

  GDouble factor = sqrt(1 + m_sign * beam_polFraction);
  complex <GDouble> zjm = 0;
  // - -> + in prod_angle
  complex <GDouble> rotateY = polar((GDouble)1., (GDouble)(-1. * prod_angle ));  

  if (m_real == 1)
	  zjm = real(amplitude * rotateY);
  if (m_real == -1) 
	  zjm = i*imag(amplitude * rotateY);

  // E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
  double kinFactor = barrierFactor(m_X, m_l, m_iso, m_ps);
  //kinFactor *= sqrt(2.*m_s + 1.) * sqrt(2.*m_l + 1.);
  factor *= kinFactor;

  userVars[uv_ampRe] = ( factor * zjm ).real();
  userVars[uv_ampIm] = ( factor * zjm ).imag();

  return;
}


/////////////////////// Amplitude Calculation //////////////////////////

complex< GDouble >
Iso_ps_refl::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  return complex< GDouble >( userVars[uv_ampRe], userVars[uv_ampIm] );
}

void Iso_ps_refl::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
Iso_ps_refl::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUIso_ps_refl_exec( dimGrid, dimBlock, GPU_AMP_ARGS );

}

#endif


