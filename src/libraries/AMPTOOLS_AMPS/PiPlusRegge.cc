
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/PiPlusRegge.h"

#include "UTILITIES/BeamProperties.h"

PiPlusRegge::PiPlusRegge( const vector< string >& args ) :
UserAmplitude< PiPlusRegge >( args )
{
	assert( args.size() == 1 );
	// Polarization plane angle (PARA = 0 and PERP = PI/2)
	PolPlane = atof( args[0].c_str() );

	// BeamProperties configuration file
	TString beamConfigFile = args[0].c_str();
	BeamProperties beamProp(beamConfigFile);
	polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
}


complex< GDouble >
PiPlusRegge::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector target  ( 0., 0., 0., 0.938);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	
	TLorentzVector cm = recoil + p1;
	TLorentzRotation cmBoost( -cm.BoostVector() );
	
	TLorentzVector beam_cm = cmBoost * beam;
	TLorentzVector target_cm = cmBoost * target;
	TLorentzVector recoil_cm = cmBoost * recoil;
	
	// phi dependence needed for polarized distribution
	TLorentzVector p1_cm = cmBoost * p1;
	GDouble phi = p1_cm.Phi() + PolPlane*TMath::Pi()/180.;
	GDouble cos2Phi = cos(2.*phi);
	
	// polarization from cobrem.F
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());
	GDouble Pgamma;
	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else Pgamma = polFrac_vs_E->GetBinContent(bin);

	GDouble t = (target - recoil).M2();
	GDouble W = exp(2.5*t);

	// hard coded beam asymmetry for all -t
	GDouble BeamSigma = 0.8;
	W *= (1 - Pgamma * BeamSigma * cos2Phi);

	return complex< GDouble > ( sqrt( fabs(W) ) );
}

