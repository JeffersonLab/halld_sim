#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/FSRootDataReaderTEM.h"
#include "TSystem.h"

using namespace std;

FSRootDataReaderTEM::FSRootDataReaderTEM( const vector< string >& args ) :
   UserDataReader< FSRootDataReaderTEM >(args),
   m_eventCounter( 0 ){


	// arguments:
	// inFileName inTreeName numParticles t_min t_max E_min E_max m_min m_max [fourMomentumPrefix]
	// inFileName inTreeName numParticles t_min t_max E_min E_max m_min m_max friendFileName friendTreeName friendBranchName [fourMomentumPrefix]
      	assert((args.size() >= 9 && args.size() <= 10) || (args.size()>=12 && args.size()<=13));
      	inFileName = args[0];
      	inTreeName = args[1];
      	m_numParticles = atoi(args[2].c_str());
      	assert (m_numParticles < 50);
	m_tMin = atof(args[3].c_str());
	m_tMax = atof(args[4].c_str());
	m_EMin = atof(args[5].c_str());
	m_EMax = atof(args[6].c_str());
	m_MMin = atof(args[7].c_str());
	m_MMax = atof(args[8].c_str());

      	TString fourMomentumPrefix = "";
      	if (args.size() == 10) fourMomentumPrefix = args[9];
      
      	TString friendFileName = "";
      	TString friendBranchName = "weight";
      	TString friendTreeName = "";
      	if (args.size() >= 12) {
        	friendFileName = args[9];
        	friendTreeName = args[10];
        	friendBranchName = args[11];
      	}
      	if (args.size() == 13) fourMomentumPrefix = args[12];

      	TH1::AddDirectory( kFALSE );
      	gSystem->Load( "libTree" );

      	ifstream fileexists( inFileName.c_str() );
      	if (fileexists){
         	m_inFile = new TFile( inFileName.c_str() );
         	m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
         	if(args.size()>=12)
            		m_inTree->AddFriend(friendTreeName, friendFileName);
      	}
      	else{
         	cout << "FSRootDataReaderTEM WARNING:  Cannot find file... " << inFileName << endl;
         	m_inFile = NULL;
         	m_inTree = NULL;
      	}

      	if(args.size()>=9)
        	cout << "Opening Tree " << inFileName << " " << inTreeName << " with " << m_numParticles << " final state particles" << endl;
      	if(args.size()==10)
        	cout << "Four-momentum prefix in FSRoot tree is " << fourMomentumPrefix << endl;
      //if(args.size()==11) //Why is this line here?
      //  cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << " " << args[3] << " " << args[4] << endl;
      	if(args.size()>=12)
        	cout << "Opening Friend Tree " << friendFileName << " " << friendBranchName << " " << friendTreeName << endl;
      	if(args.size()==13)
        	cout << "Four-momentum prefix in FSRoot tree is " << fourMomentumPrefix << endl;



      	if (m_inTree){
         	TString sEnPB = fourMomentumPrefix+"EnPB";
         	TString sPxPB = fourMomentumPrefix+"PxPB";
         	TString sPyPB = fourMomentumPrefix+"PyPB";
         	TString sPzPB = fourMomentumPrefix+"PzPB";
         	m_inTree->SetBranchAddress( sEnPB, &m_EnPB );
         	m_inTree->SetBranchAddress( sPxPB, &m_PxPB );
         	m_inTree->SetBranchAddress( sPyPB, &m_PyPB );
         	m_inTree->SetBranchAddress( sPzPB, &m_PzPB );
		
		m_numEvents = 0;
		cout << "*********************************************" << endl;
      		cout << "ROOT Data reader  -t range specified [" << m_tMin << "," << m_tMax << ")" << endl;
      		cout << "ROOT Data reader Beam E range specified [" << m_EMin << "," << m_EMax << ")" << endl;
      		cout << "ROOT Data reader  Inv. Mass range specified [" << m_MMin << "," << m_MMax << ")" << endl;
      		cout << "Total events: " <<  m_inTree->GetEntries() << endl;

         	for (unsigned int i = 0; i < m_numParticles; i++){
            		TString sI("");  sI += (i+1);
            		TString sEnPi = fourMomentumPrefix+"EnP"+sI;
            		TString sPxPi = fourMomentumPrefix+"PxP"+sI;
            		TString sPyPi = fourMomentumPrefix+"PyP"+sI;
            		TString sPzPi = fourMomentumPrefix+"PzP"+sI;
            		m_inTree->SetBranchAddress( sEnPi, &m_EnP[i] );
            		m_inTree->SetBranchAddress( sPxPi, &m_PxP[i] );
            		m_inTree->SetBranchAddress( sPyPi, &m_PyP[i] );
            		m_inTree->SetBranchAddress( sPzPi, &m_PzP[i] );
            		if(args.size()>=12)
              			m_inTree->SetBranchAddress( friendBranchName, &m_weight );
            		else
              			m_weight = 1.0;

         	}
//		while( m_eventCounter < 5 ){
		while( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
			m_inTree->GetEntry( m_eventCounter++ );
			assert( m_numParticles < Kinematics::kMaxParticles );
			
			if(checkEvent()) m_numEvents++;
		}
	
		cout << "Number of events kept    = " << m_numEvents << endl;
      		cout << "*********************************************" << endl;
      	}
}


void FSRootDataReaderTEM::resetSource(){
   	cout << "Resetting source " << inTreeName << " in " << inFileName << endl;
	m_eventCounter = 0;
}



Kinematics* FSRootDataReaderTEM::getEvent(){
	bool m_RangeSpecified = true;

	if( m_RangeSpecified == false){
		if( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
      			m_inTree->GetEntry( m_eventCounter++ );
			assert( m_numParticles < Kinematics::kMaxParticles );
     			return new Kinematics( particleList(), m_weight );	
   		}
   		else return NULL;
	}
	else{
		while( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
//		while( m_eventCounter < 5 ){
      			m_inTree->GetEntry( m_eventCounter++ );
			assert( m_numParticles < Kinematics::kMaxParticles );
 			if(checkEvent())
	    			return new Kinematics( particleList(), m_weight );	
		
		}
		return NULL;
	}
	return NULL;
}

vector<TLorentzVector> FSRootDataReaderTEM::particleList(){
	assert( m_numParticles < Kinematics::kMaxParticles );

	vector<TLorentzVector> locParticleList;
	locParticleList.push_back( TLorentzVector( m_PxPB, m_PyPB, m_PzPB, m_EnPB ) );
	for (unsigned int i = 0; i < m_numParticles; ++i){
		locParticleList.push_back( TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i] ) );
//		cout << m_PxP[i] << " " << m_PyP[i] << " " << m_PzP[i] << " " << m_EnP[i] << endl;
	}
	return locParticleList;
}

bool FSRootDataReaderTEM::checkEvent(){
	assert( m_numParticles < Kinematics::kMaxParticles );

        TLorentzVector finalstate;
        TLorentzVector recoil;

        for(unsigned int i = 0; i < m_numParticles; ++i){
//		cout << i << ": " << m_PxP[i] << " " << m_PyP[i] << " " << m_PzP[i] << " " << m_EnP[i] << endl;
           	if (i > 0 && i < m_numParticles-1){
			finalstate += TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i] );
//			cout << "Adding " << i << "th particle (with mass " << TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i]).M() << " GeV) to finalstate" << endl;
		}
            	if (i == 0 || i == m_numParticles-1){
			recoil += TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i] ); // proton (particle 0) and last particle in the list compose the recoil particle. Eventually need to modify this to allow for stable lower vertices
//			cout << "Adding " << i << "th particle (with mass " << TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i]).M() << " GeV) to recoil" << endl;
		}
	}

//	cout << "m_recoil = " << recoil.M() << endl;
//	cout << "m_finalstate = " << finalstate.M() << endl;
        // calculate -t and check if it, beam energy, and final state mass are in specified ranges
        // use the reconstructed proton
        TLorentzVector target = TLorentzVector(0, 0, 0, 0.938272);
        double tMag = fabs( (target - recoil).M2() );
        double EMag = TLorentzVector( m_PxPB, m_PyPB, m_PzPB, m_EnPB ).E();
        double MMag = finalstate.M();
        
	if(m_tMin <= tMag && tMag < m_tMax && m_EMin <= EMag && EMag < m_EMax && m_MMin <= MMag && MMag < m_MMax){
//		cout << tMag << " and " << EMag << " and " << MMag << endl;
		return true;
	}
	return false;
}


unsigned int FSRootDataReaderTEM::numEvents() const{
	if (!m_inTree) return 0;
//        return static_cast< unsigned int >( m_inTree->GetEntries() );
   	return m_numEvents;
}
