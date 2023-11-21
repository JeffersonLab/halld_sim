#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "TSystem.h"

#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

int main( int argc, char* argv[] ){

	string configfile;

	for ( int i = 1; i < argc; i++ ){
		string arg( argv[i] );

		if ( arg == "-c" ){ configfile = argv[++i]; }
	}
	
	if( configfile.size() == 0 ){
		cout << "No config file specified" << endl;
		exit(1);
	}

	ConfigFileParser parser(configfile);
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	cfgInfo->display();

	AmpToolsInterface::registerAmplitude( BreitWigner() );
	AmpToolsInterface::registerAmplitude( PhaseOffset() );
	AmpToolsInterface::registerAmplitude( OmegaDalitz() );
	AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
	AmpToolsInterface::registerAmplitude( LowerVertexDelta() );

	AmpToolsInterface ati( cfgInfo );
	vector< TLorentzVector > p4List;

	// list lab frame p4s, in the order specified in the AmpTools config file
	p4List.push_back( TLorentzVector(0, 0, 8.209470, 8.209470) ); //beam photon
	p4List.push_back( TLorentzVector(-0.291810, -0.315238, 0.695296, 1.244320) ); //recoil proton
	p4List.push_back( TLorentzVector(0.183147, 0.002392, 0.384435, 0.448137) ); //bachelor pion
	p4List.push_back( TLorentzVector(0.151104, 0.333728, 1.914300, 1.953710) ); //pi0 from omega
	p4List.push_back( TLorentzVector(-0.135341, 0.038805, 2.134320, 2.143510) ); //pi+ from omega
	p4List.push_back( TLorentzVector(0.076385, 0.004007, 3.173990, 3.177980) ); //pi- from omega
	p4List.push_back( TLorentzVector(0.016522, -0.063693, -0.092871, 0.180096) ); //recoil pion

	Kinematics kin( p4List );

	ati.printEventDetails( cfgInfo->reactionList()[0]->reactionName(), &kin );
	
	return 0;
}

