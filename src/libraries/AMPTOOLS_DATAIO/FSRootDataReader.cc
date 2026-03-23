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
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "TSystem.h"

using namespace std;

// Constructor expects one of the following argument patterns:
//
// 3 args: inFileName inTreeName numParticles
//         - Basic usage with default branch names and weight branch
//
// 4 args: inFileName inTreeName numParticles fourMomentumPrefix
//         - Adds custom prefix for four-momentum branch names
//
// 5 args: inFileName inTreeName numParticles fourMomentumPrefix weightBranchName
//         - Adds custom weight branch name (no friend tree)
//
// 6 args: inFileName inTreeName numParticles friendFileName friendTreeName weightBranchName
//         - Adds friend tree with custom weight branch name
//
// 7 args: inFileName inTreeName numParticles friendFileName friendTreeName weightBranchName fourMomentumPrefix
//         - Full specification with friend tree and custom prefix
//
FSRootDataReader::FSRootDataReader( const vector< string >& args ) :
   UserDataReader< FSRootDataReader >(args),
   m_eventCounter( 0 ){

      // Validate argument count
      assert((args.size() >= 3 && args.size() <= 7));
      
      // Parse required arguments
      string inFileName(args[0]);
      string inTreeName(args[1]);
      m_numParticles = atoi(args[2].c_str());
      assert (m_numParticles < 50);
      
      // Parse optional arguments based on count
      TString fourMomentumPrefix = "";
      TString friendFileName = "";
      TString friendTreeName = "";
      TString weightBranchName = "weight"; // default weight branch name
      
      if (args.size() == 4) {        
        fourMomentumPrefix = args[3];
      }
      else if (args.size() == 5) {        
        fourMomentumPrefix = args[3];
        weightBranchName = args[4];
      }
      else if (args.size() == 6) {
        friendFileName = args[3];
        friendTreeName = args[4];
        weightBranchName = args[5];
      }
      else if (args.size() == 7) {
        friendFileName = args[3];
        friendTreeName = args[4];
        weightBranchName = args[5];
        fourMomentumPrefix = args[6];
      }

      TH1::AddDirectory( kFALSE );
      gSystem->Load( "libTree" );

      // Open input file and tree
      ifstream fileexists( inFileName.c_str() );
      if (fileexists){
         fileexists.close();
         m_inFile = new TFile( inFileName.c_str() );
         if (!m_inFile || m_inFile->IsZombie()) {
            cout << "FSRootDataReader WARNING:  Cannot open file... " << inFileName << endl;
            m_inFile = NULL;
            m_inTree = NULL;
            return;
         }
         m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
         
         if (!m_inTree) {
            cout << "FSRootDataReader WARNING:  Cannot open tree... " << inTreeName << endl;
            m_inFile->Close();
            delete m_inFile;
            m_inFile = NULL;
            m_inTree = NULL;
            return;
         }
                  
         if(friendFileName != "" && friendTreeName != "")
            m_inTree->AddFriend(friendTreeName, friendFileName);
      }
      else{
         cout << "FSRootDataReader WARNING:  Cannot find file... " << inFileName << endl;
         m_inFile = NULL;
         m_inTree = NULL;
         return;
      }
      
      cout << "Opening Tree: " << inFileName << " " << inTreeName << " (numParticles=" << m_numParticles << ")";
      if (fourMomentumPrefix != "") cout << " fourMomentumPrefix=" << fourMomentumPrefix;
      if (friendFileName != "") cout << " friendFile=" << friendFileName << " friendTree=" << friendTreeName;
      if (weightBranchName != "weight") cout << " weightBranch=" << weightBranchName;
      cout << endl;
      if (m_inTree){
         TString sEnPB = fourMomentumPrefix+"EnPB";
         TString sPxPB = fourMomentumPrefix+"PxPB";
         TString sPyPB = fourMomentumPrefix+"PyPB";
         TString sPzPB = fourMomentumPrefix+"PzPB";
         m_inTree->SetBranchAddress( sEnPB, &m_EnPB );
         m_inTree->SetBranchAddress( sPxPB, &m_PxPB );
         m_inTree->SetBranchAddress( sPyPB, &m_PyPB );
         m_inTree->SetBranchAddress( sPzPB, &m_PzPB );
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
         }
         
         // Set up weight branch if it exists, otherwise default to 1.0
         if (friendFileName != "")
            m_inTree->SetBranchAddress( weightBranchName, &m_weight );
         else if (m_inTree->GetBranch(weightBranchName) != NULL)
            m_inTree->SetBranchAddress( weightBranchName, &m_weight );
         else
            m_weight = 1.0;
      }

   }


void FSRootDataReader::resetSource(){
   m_eventCounter = 0;
}


Kinematics* FSRootDataReader::getEvent(){
   if( m_eventCounter < numEvents() ){
      m_inTree->GetEntry( m_eventCounter++ );
      vector< TLorentzVector > particleList;
      particleList.push_back( TLorentzVector( m_PxPB, m_PyPB, m_PzPB, m_EnPB ) );
      for (unsigned int i = 0; i < m_numParticles; i++){
         particleList.push_back( TLorentzVector( m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i] ) );
      }
      return new Kinematics( particleList, m_weight );
   }
   else{
      return NULL;
   }
}


unsigned int FSRootDataReader::numEvents() const{
   if (!m_inTree) return 0;
   return static_cast< unsigned int >( m_inTree->GetEntries() );
}
