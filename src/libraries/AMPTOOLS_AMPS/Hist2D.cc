
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Hist2D.h"

Hist2D::Hist2D( const vector< string >& args ) :
UserAmplitude< Hist2D >( args )
{
	assert( args.size() == 4 );
	fileName = args[0].c_str();
	histName = args[1].c_str();
	histType = args[2].c_str();
	particleList = args[3].c_str();

	cout<<"Opening ROOT file "<<fileName.data()<<endl;
	cout<<"Model provided in histogram named "<<histName.data()<<endl;
	cout<<"Summing particle indices "<<particleList.data()<<" for invariant mass"<<endl;
	for(uint i=0; i<particleList.length(); i++) {
		string num; num += particleList[i];
		int index = atoi(num.c_str());
		cout<<index<<endl;
	}

	TFile *finput = TFile::Open(fileName.data());
	if(!finput->IsOpen()) {
		cout<<"Can't find file "<<fileName.data()<<endl;
		exit(1);
	}

	hist2D = (TH2*)finput->Get(histName.data());
	if(!hist2D) {
		cout<<"Can't find histogram "<<histName.data()<<" in file "<<fileName.data()<<endl;
		exit(1);
	}		

	// keep in memory after file is closed
        hist2D->SetDirectory(gROOT);
	finput->Close();
}


complex< GDouble >
Hist2D::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector target  ( 0., 0., 0., 0.938);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	
	// compute particle P4 sum for invariant mass
	TLorentzVector sum;
	for(uint i=0; i<particleList.length(); i++) {
		string num; num += particleList[i];
		int index = atoi(num.c_str());
		TLorentzVector particleP4 ( pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0] ); 
		sum += particleP4;
	}
	
	double beamE = beam.E();
	double t = fabs((sum - beam).M2()); 

	double userVarX = 0;
	double userVarY = sum.M();
	if(histType == "MassVsEgamma") {
		userVarX = beamE;
	}
	if(histType == "MassVst") {
		userVarX = t;
	}

	// weighted model of intensity from histogram 
	GDouble W = 0.; // initialized to zero

	int bin = hist2D->FindBin(userVarX, userVarY); // generic bin index from 2D histogram (negative value if values outside defined range)
	if(bin > 0) W = hist2D->GetBinContent(bin); 

	return complex< GDouble > ( sqrt(W) );
}
