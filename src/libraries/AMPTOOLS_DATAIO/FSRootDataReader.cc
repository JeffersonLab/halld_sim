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

#include "IUAmpTools/report.h"

const char* FSRootDataReader::kModule = "FSRootDataReader";


using namespace std;

namespace {

   // TTree::SetBranchAddress on a branch that does not exist prints a ROOT
   // error but is not fatal, which would leave the four-momentum arrays holding
   // whatever they were initialized to and let a fit run to completion on data
   // that was never read.  Collect the missing names first so that a single
   // mistake -- a wrong fourMomentumPrefix makes every branch miss at once --
   // is reported in one go rather than one branch at a time.
   //
   // Note TTree::GetBranch also searches friend trees, so this covers branches
   // supplied by an attached friend.
   void checkBranch( TTree* tree, const TString& name, vector< TString >& missing ){
      if( tree->GetBranch( name ) == NULL ) missing.push_back( name );
   }

}

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
      // a weight branch the user asked for by name is required to exist; the
      // default one is optional and falls back to a weight of 1
      bool weightBranchSpecified = false;

      if (args.size() == 4) {
        fourMomentumPrefix = args[3];
      }
      else if (args.size() == 5) {
        fourMomentumPrefix = args[3];
        weightBranchName = args[4];
        weightBranchSpecified = true;
      }
      else if (args.size() == 6) {
        friendFileName = args[3];
        friendTreeName = args[4];
        weightBranchName = args[5];
        weightBranchSpecified = true;
      }
      else if (args.size() == 7) {
        friendFileName = args[3];
        friendTreeName = args[4];
        weightBranchName = args[5];
        fourMomentumPrefix = args[6];
        weightBranchSpecified = true;
      }

      TH1::AddDirectory( kFALSE );
      gSystem->Load( "libTree" );

      // Open input file and tree
      ifstream fileexists( inFileName.c_str() );
      if (fileexists){
         fileexists.close();
         m_inFile = new TFile( inFileName.c_str() );
         if (!m_inFile || m_inFile->IsZombie()) {
            report( ERROR, kModule ) << "FSRootDataReader ERROR:  Cannot open file... " << inFileName << endl;
            m_inFile = NULL;
            m_inTree = NULL;
            return;
         }
         m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
         
         if (!m_inTree) {
            report( ERROR, kModule ) << "FSRootDataReader ERROR:  Cannot open tree... " << inTreeName << endl;
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
         report( ERROR, kModule ) << "FSRootDataReader ERROR:  Cannot find file... " << inFileName << endl;
         m_inFile = NULL;
         m_inTree = NULL;
         return;
      }
      
     report( DEBUG, kModule ) << "Opening Tree: " << inFileName << " " << inTreeName << " (numParticles=" << m_numParticles << ")";
     if (fourMomentumPrefix != "") report( DEBUG, kModule ) << " fourMomentumPrefix=" << fourMomentumPrefix;
     if (friendFileName != "") report( DEBUG, kModule ) << " friendFile=" << friendFileName << " friendTree=" << friendTreeName;
     if (weightBranchName != "weight") report( DEBUG, kModule ) << " weightBranch=" << weightBranchName;
     report( DEBUG, kModule ) << endl;

      if (m_inTree){
         TString sEnPB = fourMomentumPrefix+"EnPB";
         TString sPxPB = fourMomentumPrefix+"PxPB";
         TString sPyPB = fourMomentumPrefix+"PyPB";
         TString sPzPB = fourMomentumPrefix+"PzPB";

         // check that every branch we are about to read actually exists before
         // binding any of them -- see the note on checkBranch above
         vector< TString > missing;
         checkBranch( m_inTree, sEnPB, missing );
         checkBranch( m_inTree, sPxPB, missing );
         checkBranch( m_inTree, sPyPB, missing );
         checkBranch( m_inTree, sPzPB, missing );
         for (unsigned int i = 0; i < m_numParticles; i++){
            TString sI("");  sI += (i+1);
            checkBranch( m_inTree, fourMomentumPrefix+"EnP"+sI, missing );
            checkBranch( m_inTree, fourMomentumPrefix+"PxP"+sI, missing );
            checkBranch( m_inTree, fourMomentumPrefix+"PyP"+sI, missing );
            checkBranch( m_inTree, fourMomentumPrefix+"PzP"+sI, missing );
         }
         if (weightBranchSpecified) checkBranch( m_inTree, weightBranchName, missing );

         if (!missing.empty()){
            report( ERROR, kModule ) << "FSRootDataReader ERROR:  " << missing.size()
               << " branch(es) not found in tree " << inTreeName
               << " of file " << inFileName << ":" << endl;
            for (unsigned int i = 0; i < missing.size(); i++)
               report( ERROR, kModule ) << "     " << missing[i] << endl;
            report( ERROR, kModule ) << "  check the arguments:  numParticles = "
               << m_numParticles << ", fourMomentumPrefix = \"" << fourMomentumPrefix
               << "\", weightBranchName = \"" << weightBranchName << "\"" << endl;
            exit( 1 );
         }

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

         // a weight branch named by the user is bound above having been checked;
         // the default one is used only if the tree happens to provide it
         if (weightBranchSpecified || m_inTree->GetBranch(weightBranchName) != NULL)
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
