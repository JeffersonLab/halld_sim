
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/SinglePS.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/decayAngles.h"

#include "UTILITIES/BeamProperties.h"

SinglePS::SinglePS( const vector< string >& args ) :
UserAmplitude< SinglePS >( args )
{
  m_r = atoi( args[0].c_str() ); // real (+1) or imaginary (-1)
  m_s = atoi( args[1].c_str() ); // sign for polarization in amplitude

  // default polarization information stored in tree
  m_polInTree = true;

  // 2 possibilities to initialize this amplitude:
  // (with <r>: +1/-1 for real/imaginary part; <s>: +1/-1 sign in P_gamma term)

  // loop over any additional amplitude arguments to change defaults
  for( uint ioption = 2; ioption < args.size(); ioption++ ) {
	  TString option = args[ioption].c_str();

	  // polarization provided in configuration file
	  if( ioption == 2 && option.IsFloat() ) {
		  m_polInTree = false;
		  polAngle = atof( args[2].c_str() );
	  
		  TString polOption = args[3].c_str();
		  if( polOption.IsFloat() ) polFraction = atof( polOption.Data() );
		  else if(polOption.Contains(".root")) {
			  polFraction = 0.;
			  TFile* f = new TFile( polOption );
			  polFrac_vs_E = (TH1D*)f->Get( args[4].c_str() );
			  assert( polFrac_vs_E != NULL );
		  }
		  else {
			  cout << "ERROR: SinglePS beam polarization not set" <<endl;
			  assert(0);
		  }
	  }
  }

  // make sure values are reasonable
  // m_r = +1 for real
  // m_r = -1 for imag
  assert( abs( m_r ) == 1 );
  // m_s = +1 for 1 + Pgamma
  // m_s = -1 for 1 - Pgamma
  assert( abs( m_s ) == 1 );
  
}

void
SinglePS::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector beam;
  TVector3 eps;
  double beam_polFraction;
  double beam_polAngle;

  if(m_polInTree){
    beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
    eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // beam polarization vector;

    beam_polFraction = eps.Mag();
    beam_polAngle = eps.Phi()*TMath::RadToDeg();
  }
  else {
    beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
    beam_polAngle = polAngle;
    
    if(polFraction > 0.) { // for fitting with fixed polarization
	    beam_polFraction = polFraction;
    }
    else { // for fitting with polarization vs E_gamma from input histogram 
	    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
	    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		    beam_polFraction = 0.;
	    } else 
		    beam_polFraction = polFrac_vs_E->GetBinContent(bin);
    }
  }
  
  TLorentzVector ps ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); // pi-

  TLorentzVector ps_recoil(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]); // pi+
  TLorentzVector proton_recoil(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]); // proton

  TLorentzVector recoil = ps_recoil + proton_recoil;

  //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  TLorentzVector target(0,0,0,0.938);
  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;

  // set beam polarization angle to 0 degrees; apply diamond orientation in calcAmplitude
  double phiProd = getPhiProd( 0., recoil, beam, target, 2, false ); 

  userVars[uv_prod_Phi] = phiProd;

  userVars[uv_beam_polFraction] = beam_polFraction;
  userVars[uv_beam_polAngle] = beam_polAngle;

  return;
}


////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
SinglePS::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{
  GDouble prod_angle = userVars[uv_prod_Phi];
  GDouble beam_polFraction = userVars[uv_beam_polFraction];
  GDouble beam_polAngle = userVars[uv_beam_polAngle];

  complex <GDouble> amplitude(0,0);
  complex <GDouble> i(0,1);

  GDouble Factor = sqrt(1 + m_s * beam_polFraction);
  
  complex< GDouble > rotateY = polar( (GDouble)1., (GDouble)(-1.*(prod_angle + beam_polAngle*TMath::DegToRad())) ); // - -> + in prod_angle and polAngle summing
  
  if( m_r == 1 )
	  amplitude = real( rotateY );
  if( m_r == -1 ) 
	  amplitude = i*imag( rotateY );

  // E852 Nozar thesis has sqrt(2*s+1)*sqrt(2*l+1)*F_l(p_omega)*sqrt(omega)
//  double kinFactor = barrierFactor(MX, m_l, MVec, MPs);
  //kinFactor *= sqrt(3.) * sqrt(2.*m_l + 1.);
//  Factor *= kinFactor;

  return complex< GDouble >( static_cast< GDouble>( Factor ) * amplitude );
}


void SinglePS::updatePar( const AmpParameter& par ){

  // could do expensive calculations here on parameter updates  
}


#ifdef GPU_ACCELERATION

void
SinglePS::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUSinglePS_exec( dimGrid, dimBlock, GPU_AMP_ARGS, m_r, m_s );

}

#endif


