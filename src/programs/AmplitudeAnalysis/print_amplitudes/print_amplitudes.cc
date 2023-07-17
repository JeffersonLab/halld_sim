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
	AmpToolsInterface::registerAmplitude( OmegaDalitz() );
	AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
	AmpToolsInterface::registerAmplitude( LowerVertexDelta() );

	AmpToolsInterface ati( cfgInfo );
	vector< TLorentzVector > p4List;
	
	// input Vincent's events
	p4List.push_back( TLorentzVector(0, 0, 1.9437, 1.9437) ); //beam photon
	p4List.push_back( TLorentzVector(0.3300, -0.5918, -1.7263, 2.0782) ); //recoil proton
	p4List.push_back( TLorentzVector(-0.0548, 0.0983, 0.2882, 0.3375) ); //bachelor pion
	p4List.push_back( TLorentzVector(-0.0371, 0.1172, 0.3739, 0.4161) ); //pi0 from omega
	p4List.push_back( TLorentzVector(-0.2106, 0.0912, 0.2905, 0.3941) ); //pi+ from omega
	p4List.push_back( TLorentzVector(0.0291, 0.1849, 0.4849, 0.5370) ); //pi- from omega
	p4List.push_back( TLorentzVector(-0.0566, 0.1003, 0.2888, 0.3390) ); //recoil pion
	
	Kinematics kin( p4List );

//	ati.printEventDetails( omegapi, &kin );
	ati.printEventDetails( cfgInfo->reactionList()[0]->reactionName(), &kin );
	
	return 0;
}

