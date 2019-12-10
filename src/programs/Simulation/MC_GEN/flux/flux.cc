 /********************************************************
 *  Usage: flux2ascii  < input.gen
 ********************************************************
 *                      * grab beam config files        
 *  flux.cc             * invoke beamProperties class        
 *                      * translate the beam energy distribution into ASCII file
 ********************************************************
 * 
 * created by:  Hao Li
 *              Carnegie Mellon University
 *              06-Dec-2019
 ******************************************************** */
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string>
#include "UTILITIES/BeamProperties.h"

#include <particleType.h>
//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
	TString beamConfigFilename;
	TH1D *cobrem_vs_E = NULL;
	// read in beamConfigFilename somehow
	if (beamConfigFilename.Contains(".conf")) 
	{
		BeamProperties beamProp(beamConfigFilename.Data());
		cobrem_vs_E = (TH1D*)beamProp.GetFlux();
	}
	if(beamConfigFile.Length() == 0) 
	{
		cout<<"WARNING: Couldn't find beam configuration file -- write local
		version"<<endl;
		beamConfigFile = "local_beam.conf";
		ofstream locBeamConfigFile;
		locBeamConfigFile.open(beamConfigFile.Data());
		locBeamConfigFile<<"ElectronBeamEnergy "<<beamMaxE<<endl; // electron
		beam energy
		locBeamConfigFile<<"CoherentPeakEnergy "<<beamPeakE<<endl; // coherent
		peak energy
		locBeamConfigFile<<"PhotonBeamLowEnergy "<<beamLowE<<endl; // photon
		beam low energy
		locBeamConfigFile<<"PhotonBeamHighEnergy "<<beamHighE<<endl; // photon
		beam high energy
		locBeamConfigFile.close();
	}

	TString outName = "beamProfile.ascii";
	ofstream outfile;
	outfile.open(outName, ios::out);

	for(int i = 1; i< cobrem_vs_E->GetNbinsX()+1;i++)
	{
		outfile << cobrem_vs_E->GetBinCenter(i) << " " << cobrem_vs_E->GetBinContent(i) << endl;
	}
	outfile.close();
}