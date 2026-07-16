//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Here are some useful notes to understand how this amplitude works and how to
// use it in your AmpTools configuration file:
//
// PURPOSE:
//   Calculates the lineary polarized photoproduction amplitude for a system 
//   composed of a vector (V) and pseudoscalar (ps)
//
// PHYSICS SUMMARY:
//   The amplitude uses the helicity formalism (see [GlueX-doc-4858]): 
//   a resonance X with total spin J, spin projection M, and partial wave L 
//   decays into V and ps. The angular distribution is described by 
//   Wigner D-functions and Clebsch-Gordan coefficients. A barrier factor
//   accounts for the centrifugal suppression near threshold.
//
//   Beam polarization information is required in two places:
//   - The fraction linear polarization.
//   - The angle between the polarization angle and the production plane
//
//   The default behavior is that the vector decays into two pseudoscalars
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CONSTRUCTOR ARGUMENTS (for the AmpTools config file):
// We support two formats for optional arguments: a positional input format 
// and a key=value format
// 
// REQUIRED--the first 5 arguments are required and positional:
//   [0] m_j    - Total spin (J) of the resonance X (integer >= 0)
//   [1] m_m    - Spin projection (m) of X; must satisfy |m| <= J
//   [2] m_l    - Orbital angular momentum (l) between V and ps (integer >= 0)
//   [3] m_r    - +1 (r)eal part or -1 imaginary part of the amplitude
//   [4] m_s    - (S)ign of the polarization fraction +1 for (1 + P_gamma) or
//                 -1 for (1 - P_gamma)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OPTIONAL--positional arguments:
//   [5] polAngle     Beam polarization angle. Degrees will be converted 
//                    to radians
//   [6] polFraction  Fixed polarization fraction [0, 1],
//                    OR a path to a .root file containing a TH1D histogram
//                    of polarization vs beam energy
//   [7] histName     (required if [6] is a .root file) - Name of the TH1D
//                    inside the root file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OPTIONAL--key=value arguments:
//   polAngle=<val>       Beam polarization angle
//   polFraction=<val>    Fixed polarization fraction [0, 1]
//   polFile=<path.root>  path to a .root file containing a TH1D histogram
//                        of polarization vs beam energy
//   polHist=<histName>   Must be used if polFile is used
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Note: If NO Optional arguments are given, the code assumes that you have
// stored the polarization information in the Px and Py components of the 
// beam. This is useful since the acceptance within a run period is assumed
// to be constant and reduces the data set by a factor of 4.
// When constructing your tree save the beam 4-vector as follows:
// Px = polFraction*cos(polAngle)  *Make sure the angle is in radians
// Py = polFraction*sin(polAngle)  *Make sure the angle is in radians
// Pz = Regular beam Pz
// E  = Regular beam energy        *Note that Pz = E
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The following are keywords to turn-on features on the amplitude. They must 
// be placed at the end of optional arguments if using the positional format, 
// or can be placed anywhere if using the key=value format:
//   [X] "omega3pi"        V decays will be described as omega -> pi+ pi- pi0
//   [Y] "omegagpi0"       V decays will be described as omega -> gamma pi0
//   [Y+1] gHelicity       When omegagpi0 is used, one must specify the 
//                         helicity of the bachelor photon (+1 or -1)
//   gHelicity=<+1|-1>     Alternative key=value format for the photon helicity
///  "nobarrier"           Turn off the barrier factors
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Examples of usage in config file:
//  For a 1S m=-1 wave from an omega->3pi that belongs to the imaginary part of
//  the amplitude sum and has a negative sign in front of the polarization term,
//  you would write:
//  If polarization info stored in the beam P4:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 1 -1 0 -1 -1 omega3pi
//  If polarization angle and fraction are fixed:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 1 -1 0 -1 -1 90 .35 omega3pi
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 1 -1 0 -1 -1 polFraction=0.35 polAngle=90 omega3pi
//  If polarization vs E_gamma from histogram:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 1 -1 0 -1 -1 90 pathToFile.root histName omega3pi
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 1 -1 0 -1 -1 polAngle=90 polFile=pathToFile.root polHist=histName omega3pi
//  For a 2P m=0 wave from an omega->gpi that belongs to the real part of
//  the amplitude sum,has a negative sign in front of the polarization term, and
//  a photon helicity -1, you would write:
//  If polarization info stored in the beam P4:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 omegagpi -1
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 omegagpi gHelicity=-1
//  If polarization angle and fraction are fixed:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 90 .35 omegagpi -1
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 polFraction=0.35 polAngle=90 omegagpi gHelicity=-1
//  If polarization vs E_gamma from histogram:
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 90 pathToFile.root histName omegagpi -1
//  amplitude reactionName::SumName::WaveName Vec_ps_refl 2 0 1 1 -1 polAngle=90 polFile=pathToFile.root polHist=histName omegagpi gHelicity=-1
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// When the reaction is specified in the config file, the following format is 
// assumed: 
// reaction reactionName Beam Proton Ps VecDaught1 VecDaught2 [VecDaught3]
//
// This order should be reflected in the way the particles are saved in the 
// input file. Mistakes deviating from this assumption will not trigger
// an error and will lead to incorrect results.
// Assumed ordering of particles in the array from the input file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/vecPsAngles.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "IUAmpTools/report.h"

#include "UTILITIES/BeamProperties.h"

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//....oooOO0OOooo........ Helper Functions ........oooOO0OOooo.....
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	Vec_ps_refl::VecPsReflArgs Vec_ps_refl::parsedArgs(
                    const std::vector<std::string>& args, 
									  const std::string& context){
	  VecPsReflArgs inputArgs;
	
	  // Required positional args 
	  inputArgs.j = static_cast<int>( parseValidatedNumber( "J",            args[0], context ) );
	  inputArgs.m = static_cast<int>( parseValidatedNumber( "M",            args[1], context ) );
	  inputArgs.l = static_cast<int>( parseValidatedNumber( "L",            args[2], context ) );
	  inputArgs.r = static_cast<int>( parseValidatedNumber( "Re/Im",        args[3], context ) );
	  inputArgs.s = static_cast<int>( parseValidatedNumber( "P_gamma sign", args[4], context ) );
	
	  // Physics validation of required args
	  if( abs(inputArgs.m) > inputArgs.j )
		throw std::invalid_argument(
		  "[ " + context + " ] |M| must be <= J, got J=" + args[0] + " M=" + args[1] );
	  if( abs(inputArgs.r) != 1 )
		throw std::invalid_argument(
		  "[ " + context + " ] Re/Im flag must be +1 or -1, got " + args[3] );
	  if( abs(inputArgs.s) != 1 )
		throw std::invalid_argument(
		  "[ " + context + " ] P_gamma sign must be +1 or -1, got " + args[4] );
	
	  // If no optional args, return now
	  if( args.size() <= 5 ) return inputArgs;
	
	  // Detect format by checking if args[5] contains '='
	  bool isKeyValue = std::any_of( args.begin() + 5, args.end(),
		  []( const std::string& s ){ return s.find('=') != std::string::npos; } );
	
	  // Lambda function to detect bare flags (valid in both formats)
	  auto isBareFlag = []( const std::string& s ){
		return s == "omega3pi" || s == "omegagpi0" || s == "nobarrier";
	  };
	
	  if( !isKeyValue ){
		size_t i = 5;
	
		// arg[5]: polarization angle
		if( i < args.size() && !isBareFlag(args[i]) && args[i].find('=') == std::string::npos ){
		  inputArgs.polAngle  = parseValidatedNumber( "polAngle", args[i++], context );
		  inputArgs.polInfoInPhotonP4 = false;
		}
	
		// arg[6]: fixed fraction or .root file path
		if( i < args.size() && !isBareFlag(args[i]) && args[i].find('=') == std::string::npos ){
		  if( args[i].find(".root") != std::string::npos ){
			inputArgs.polFile = args[i++];
			if( i >= args.size() || isBareFlag(args[i]) )
			  throw std::invalid_argument(
				"[ " + context + " ] .root polarization file requires a histogram name immediately after it" );
			inputArgs.polHist = args[i++];
		  } else {
			inputArgs.polFraction = parseValidatedNumber( "polFraction", args[i++], context );
		  }
		}
	
		// Remaining: bare flags, and gHelicity bare integer after omegagpi0
		while( i < args.size() ){
		  if( args[i] == "omega3pi" ){
			inputArgs.omega3pi = true;
			i++;
		  }
		  else if( args[i] == "omegagpi0" ){
			inputArgs.omegagpi0 = true;
			if( i + 1 >= args.size() )
			  throw std::invalid_argument(
				"[ " + context + " ] omegagpi0 requires a photon helicity value (+1 or -1) immediately after it" );
			inputArgs.gHelicity = static_cast<int>( parseValidatedNumber( "gHelicity", args[i+1], context ) );
			i += 2;
		  }
		  else if( args[i] == "nobarrier" ){
			inputArgs.noBarrier = true;
			i++;
		  }
		  else{
			throw std::invalid_argument(
			  "[ " + context + " ] unrecognized positional argument '" + args[i] + "'" );
		  }
		}
	  } 
	  else { // key=value format
		for( size_t i = 5; i < args.size(); i++ ){
		  const std::string& arg = args[i];
	
		  if( arg == "omega3pi"  ){ inputArgs.omega3pi  = true; continue; }
		  if( arg == "omegagpi0" ){ inputArgs.omegagpi0 = true; continue; }
		  if( arg == "nobarrier" ){ inputArgs.noBarrier = true; continue; }
		  auto sep = arg.find('='); //  '=' is used as separator
		  if( sep == std::string::npos )
			throw std::invalid_argument(
			  "[ " + context + " ] unrecognized argument '" + arg +
			  "' (expected key=value pair or known flag)" );
	
		  const std::string key   = arg.substr(0, sep);
		  const std::string value = arg.substr(sep + 1);
	
		  if( key == "polAngle" ){
			inputArgs.polAngle  = parseValidatedNumber( "polAngle",    value, context );
			inputArgs.polInfoInPhotonP4 = false;
		  }
		  else if( key == "polFraction" ){
			inputArgs.polFraction = parseValidatedNumber( "polFraction", value, context );
		  }
		  else if( key == "polFile" ){ inputArgs.polFile = value; }
		  else if( key == "polHist" ){ inputArgs.polHist = value; }
		  else if( key == "gHelicity" ){
			inputArgs.gHelicity = static_cast<int>( parseValidatedNumber( "gHelicity", value, context ) );
		  }
		  else{
			throw std::invalid_argument(
			  "[ " + context + " ] unrecognized key '" + key + "'" );
		  }
		}
	  }
	
	  // Post-parse validation
	  if( inputArgs.polFile != "" && inputArgs.polHist == "" )
		throw std::invalid_argument(
		  "[ " + context + " ] polFile='" + inputArgs.polFile +
		  "' requires polHist=<histName> to also be specified" );
	
	  if( inputArgs.polHist != "" && inputArgs.polFile == "" )
		throw std::invalid_argument(
		  "[ " + context + " ] polHist='" + inputArgs.polHist +
		  "' was specified without polFile=<path.root>" );
	
	  if( inputArgs.polFraction >= 0.0 && inputArgs.polFile != "" )
		throw std::invalid_argument(
		  "[ " + context + " ] polFraction and polFile are mutually exclusive; specify only one" );
		  
	  if( inputArgs.omegagpi0 && abs(inputArgs.gHelicity) != 1 )
		throw std::invalid_argument(
		  "[ " + context + " ] omegagpi0 requires gHelicity=+1 or gHelicity=-1" );
		  
	  if( inputArgs.gHelicity != 0 && !inputArgs.omegagpi0 )
		throw std::invalid_argument(
		  "[ " + context + " ] gHelicity was specified but omegagpi0 flag is missing" );
	
	  if( inputArgs.omega3pi && inputArgs.omegagpi0 )
		throw std::invalid_argument(
		  "[ " + context + " ] omega3pi and omegagpi0 are mutually exclusive; specify only one" );
	
	  if( !inputArgs.polInfoInPhotonP4 && inputArgs.polFraction < 0.0 && inputArgs.polFile == "" )
		throw std::invalid_argument(
		  "[ " + context + " ] polAngle was given but neither polFraction nor polFile+polHist were provided" );
	
	  return inputArgs;
	}

Vec_ps_refl::Vec_ps_refl( const vector< string >& args ) :
UserAmplitude< Vec_ps_refl >( args ){
  VecPsReflArgs inputArgs = parsedArgs( args );
  m_j                  = inputArgs.j;
  m_m                  = inputArgs.m;
  m_l                  = inputArgs.l;
  m_r                  = inputArgs.r;
  m_s                  = inputArgs.s;
  
  m_polInfoInPhotonP4  = inputArgs.polInfoInPhotonP4;
  m_polAngle           = inputArgs.polAngle;
  m_polFraction        = inputArgs.polFraction;
  m_polFracVsE         = nullptr;

  m_3pi                = inputArgs.omega3pi;
  m_gpi0               = inputArgs.omegagpi0;
  m_ghel               = inputArgs.gHelicity;
  m_noBarrier          = inputArgs.noBarrier;
  

  // Open polarization file and retrieve histogram
  if( inputArgs.polFile != "" ){
    TFile* f = new TFile( inputArgs.polFile.c_str() );
    if( !f || f->IsZombie() )
      throw std::runtime_error(
        "[ Vec_ps_refl ]: could not open polarization file '" + inputArgs.polFile + "'" );
    m_polFracVsE = (TH1D*)f->Get( inputArgs.polHist.c_str() );
    if( !m_polFracVsE )
      throw std::runtime_error(
        "[ Vec_ps_refl ]: histogram '" + inputArgs.polHist +
        "' not found in '" + inputArgs.polFile + "'" );
    // Detach histogram from file so it persists after file is closed
    m_polFracVsE->SetDirectory(0);
    f->Close();
    delete f;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//....oooOO0OOooo........ Main Functions ........oooOO0OOooo.....
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void
Vec_ps_refl::calcUserVars( GDouble** pKin, GDouble* userVars ) const{

  TLorentzVector beam;
  TVector3 eps;              // beam polarization vector (eps)ilon
  GDouble beamPolFraction;
  GDouble beamPolAngle;

  if(m_polInfoInPhotonP4){
    // When pol info is stored in the photon 4-vector in the tree
    // the energy (pKin[0][0]) and the pz (pKin[0][3]) are used as normal 
    beam.SetPxPyPzE( 0.0, 0.0, pKin[0][3], pKin[0][0]);
    // while the px and py components store the pol info
    // The values should be stored as px = polFraction*cos(polAngle) 
    // and py = polFraction*sin(polAngle)
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.0); 
    beamPolFraction = eps.Mag();
    beamPolAngle = eps.Phi();
  }
  else{
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    beamPolAngle = m_polAngle;
    // for fixed polarization fraction
    if(m_polFracVsE == nullptr) 
	    beamPolFraction = m_polFraction;
    else{ // for fitting with polarization vs E_gamma from input histogram 
	    int bin = m_polFracVsE->GetXaxis()->FindBin(pKin[0][0]);
	    if (bin == 0 || bin > m_polFracVsE->GetXaxis()->GetNbins()){
		    throw std::runtime_error(
        "[ Vec_ps_refl ]: energy " + std::to_string(pKin[0][0]) + 
        " outside histogram range" );
	    } else
		    beamPolFraction = m_polFracVsE->GetBinContent(bin);
    }
  }
  
  TLorentzVector recoil( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 

  // Fill in four-vectors for final state particles
  // 1st after proton is always the pseudoscalar meson
  TLorentzVector ps(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); 
  // Compute vector meson from its decay products
  // Make sure the order of daughters is correct in the config file!
  TLorentzVector vec, vecDaught1, vecDaught2; 

  if(m_3pi){
    // Omega ps proton, omega -> 3pi (6 particles):
    // beam proton ps pi0 pip pim
    // Omega pi- Delta++, omega -> 3pi (6 particles):
    // beam delta ps pi0 pip pim
	  TLorentzVector pi0(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  TLorentzVector pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  TLorentzVector pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
	  vec = pi0 + pip + pim;
	  vecDaught1 = pip;
	  vecDaught2 = pim;
  }
  else{
    // Omega ps proton, omega -> gpi0 (5 particles):
    // beam proton ps pi0 gamma
    // Omega pi- Delta++, omega -> gpi (5 particles):
    // beam delta ps pi0 gamma
    // Vec(KK) ps proton, (5 particles)
    // beam proton ps K K
	// Vec(pipi) pi- Delta++, (5 particles)
    // beam proton ps pi pi
	// Vec(Kpi) K+ Lambda, (5 particles)
    // beam proton ps K pi
      vecDaught1 = TLorentzVector(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
	  vecDaught2 = TLorentzVector(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
	  vec = vecDaught1 + vecDaught2;
  }

  // Final meson system P4
  TLorentzVector xMeson = vec + ps;

  ///////////////// Boost Particles and Get Angles/////////////////////

  TLorentzVector target(0,0,0,0.938); // proton at rest
  TLorentzVector beamTarget = beam + target;

  // Calculate decay angles for X in helicity frame (same for all vectors)
  // Change getXDecayAngles to get Gottfried-Jackson angles if needed
  // Note: it also calculates the production angle
  vector <double> xDecayAngles = getXDecayAngles( beamPolAngle, beam, beamTarget, xMeson, vec);

  // Calculate vector decay angles (unique for each vector)
  vector <double> vectorDecayAngles;
  if(m_3pi){
    vectorDecayAngles = getVectorDecayAngles( beamTarget, xMeson, vec,
                                              vecDaught1, vecDaught2);
  }
  else{
    vectorDecayAngles = getVectorDecayAngles( beamTarget, xMeson, vec,
                                        vecDaught1, TLorentzVector(0,0,0,0));
  }

  GDouble cosTheta = TMath::Cos(xDecayAngles[0]);
  GDouble phi = xDecayAngles[1];
  GDouble prodAngle = xDecayAngles[2]; // bigPhi
  GDouble cosThetaH = TMath::Cos(vectorDecayAngles[0]);
  GDouble phiH = vectorDecayAngles[1];
  GDouble xMesonMass = xMeson.M();
  GDouble vecMass = vec.M();
  GDouble psMass = ps.M();

  complex <GDouble> amplitude(0,0);
  static const complex <GDouble> i(0,1);

  if(m_gpi0){ // radiative omega decay requires a handling of the photon helicity
    for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
      GDouble helAmp = clebschGordan(m_l, 1, 0, lambda, m_j, lambda);
      if(lambda==0){
	      amplitude += conj(wignerD(m_j, m_m, lambda, cosTheta, phi)) *
	                   helAmp * conj(wignerD(1, lambda, m_ghel, cosThetaH, phiH)) *
	                   (m_ghel*1.0); // m_ghel is int
      }
      else{
	      amplitude += conj(wignerD(m_j, m_m, lambda, cosTheta, phi)) *
	                    helAmp * conj(wignerD(1, -1*lambda, m_ghel, cosThetaH, phiH)) *
                      -1.0;
	      // the power of lambda is irrelevant for lambda =/= 0: (-1)^lambda = -1
      }
    }
  }
  else{ // for any other vector decay
    for (int lambda = -1; lambda <= 1; lambda++) { // sum over vector helicity
      GDouble helAmp = clebschGordan(m_l, 1, 0, lambda, m_j, lambda);
            amplitude += conj(wignerD(m_j, m_m, lambda, cosTheta, phi)) *
                        helAmp * conj(wignerD(1, lambda, 0, cosThetaH, phiH));
    }
  }

  // The amplitude is multiplied by a factor, either sqrt(1 + -P_gamma) or
  // sqrt(1 + P_gamma) depending on the what sum is being calculated
  GDouble factor = sqrt(1 + m_s * beamPolFraction);
  // The result of the function that depends on the angles is stored in zjm
  complex <GDouble> zjm = 0;
  // A - sign translates to + in the prod_angle
  // This is because reflectivity conventions introduce exp(-iPhi)
  complex <GDouble> rotateY = polar((GDouble)1., (GDouble)(-1. * prodAngle ));  

  if (m_r == 1)
	  zjm = real(amplitude * rotateY);
  if (m_r == -1) 
	  zjm = i*imag(amplitude * rotateY);

  if( !m_noBarrier ){
  GDouble kinFactor = barrierFactor(xMesonMass, m_l, vecMass, psMass);
  factor *= kinFactor;
  // Alternatives:
  // E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
  // kinFactor *= sqrt(3.) * sqrt(2.*m_l + 1.);
  }
  
  userVars[uv_ampRe] = ( factor * zjm ).real();
  userVars[uv_ampIm] = ( factor * zjm ).imag();

  return;
}

/////////////////////// Amplitude Calculation //////////////////////////

complex< GDouble >
Vec_ps_refl::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  return complex< GDouble >( userVars[uv_ampRe], userVars[uv_ampIm] );
}

void Vec_ps_refl::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
Vec_ps_refl::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUVec_ps_refl_exec( dimGrid, dimBlock, GPU_AMP_ARGS );

}

#endif


