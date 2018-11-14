
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Compton.h"
#include "UTILITIES/BeamProperties.h"

Compton::Compton( const vector< string >& args ) :
UserAmplitude< Compton >( args )
{
	assert( args.size() == 1 );

	// BeamProperties configuration file
        TString beamConfigFile = args[0].c_str();
        BeamProperties beamProp(beamConfigFile);
        polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
        polAngle = beamProp.GetPolAngle();
}


complex< GDouble >
Compton::calcAmplitude( GDouble** pKin ) const {
  
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
	GDouble phi = p1_cm.Phi() + polAngle*TMath::Pi()/180.;
	GDouble cos2Phi = cos(2.*phi);
	
	// polarization from cobrem.F
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());
	GDouble Pgamma;
	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else Pgamma = polFrac_vs_E->GetBinContent(bin);	

	// factors needed to calculate cross section from model
	GDouble s = cm.M2();
	GDouble t = (recoil - target).M2();
	if(fabs(t) < 0.15) return 0.;

	// model parameters from PAC42 proposal PR12-14-003
	GDouble s_0 = 10.92;
	GDouble t_0 = 2.61;
	GDouble W = 0.0702 * pow(s_0/s, 2) * pow(t_0/t, 4);

	// include polarization effects later
	GDouble BeamSigma = 0.1;	
	W *= (1 - Pgamma * BeamSigma * cos2Phi);

	return complex< GDouble > ( sqrt( fabs(W) ) );
}

