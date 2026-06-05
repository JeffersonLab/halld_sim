#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <set>  // Added for multiset
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/FSRootDataReaderBootstrap.h"
#include "TSystem.h"

#include "IUAmpTools/report.h"

const char* FSRootDataReaderBootstrap::kModule = "FSRootDataReaderBootstrap";

// Constructor expects one of the following argument patterns:
//
// 3 args: inFileName inTreeName numParticles
//         - Basic usage with default branch names, weight branch, and seed=0
//
// 4 args: inFileName inTreeName numParticles {fourMomentumPrefix|randSeed}
//         - If arg[3] is integer: used as randSeed
//         - Otherwise: used as fourMomentumPrefix
//
// 5 args: inFileName inTreeName numParticles fourMomentumPrefix randSeed
//         - Adds custom prefix and random seed
//
// 6 args: inFileName inTreeName numParticles friendFileName friendTreeName weightBranchName
//
// 7 args: inFileName inTreeName numParticles friendFileName friendTreeName weightBranchName {randSeed|fourMomentumPrefix}
//         - If args[6] is an integer, it's interpreted as a random seed
//         - otherwise it's a four momentum prefix
/
// 8 args: inFileName inTreeName numParticles friendFileName friendTreeName weightBranchName fourMomentumPrefix randSeed
//         - Full specification with friend tree, custom prefix, and random seed
//
FSRootDataReaderBootstrap::FSRootDataReaderBootstrap( const vector< string >& args ) :
   UserDataReader< FSRootDataReaderBootstrap >(args),
   m_eventCounter( 0 ){

      // Validate argument count
      assert((args.size() >= 3 && args.size() <= 8));
      
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
      int randSeed = 0;

      auto isInteger = [](const string& s) {
          return !s.empty() && s.find_first_not_of("-0123456789") == string::npos;
      };

      if (args.size() == 4) {
          if (isInteger(args[3])) {
              randSeed = atoi(args[3].c_str());
          } else {
              fourMomentumPrefix = args[3];
          }
      } 
      else if (args.size() == 5) {
          fourMomentumPrefix = args[3];
          randSeed = atoi(args[4].c_str());
      } 
      else if (args.size() >= 6) {
          friendFileName = args[3];
          friendTreeName = args[4];
          weightBranchName = args[5];
        
        if (args.size() == 7) {
          if (isInteger(args[6])) {
            randSeed = atoi(args[6].c_str());
          } else {
            fourMomentumPrefix = args[6];
          }
        } else if (args.size() == 8) {
          fourMomentumPrefix = args[6];
          randSeed = atoi(args[7].c_str());
        }
      }

      // Initialize random number generator with provided seed
      m_randGenerator = new TRandom3(randSeed);

      // Compact printing of parsed arguments
     report( DEBUG, kModule ) << "FSRootDataReaderBootstrap initialized with:\n"
          << "  inFileName:        " << inFileName << "\n"
          << "  inTreeName:        " << inTreeName << "\n"
          << "  numParticles:      " << m_numParticles << "\n"
          << "  friendFileName:    " << (friendFileName.Length() == 0 ? "N/A" : friendFileName.Data()) << "\n"
          << "  friendTreeName:    " << (friendTreeName.Length() == 0 ? "N/A" : friendTreeName.Data()) << "\n"
          << "  weightBranchName:  " << (weightBranchName.Length() == 0 ? "N/A" : weightBranchName.Data()) << "\n"
          << "  fourMomentumPrefix:" << (fourMomentumPrefix.Length() == 0 ? "N/A" : fourMomentumPrefix.Data()) << "\n"
          << "  randSeed:          " << randSeed << "\n"
          << std::endl;


     report( NOTICE, kModule ) << "******************** NOTICE ************************" << endl;
     report( NOTICE, kModule ) << "*  You are using the boostrap data reader, which   *" << endl;
     report( NOTICE, kModule ) << "*  should only be used for evaluating errors.      *" << endl;
     report( NOTICE, kModule ) << "*  The results with different seeds will be random *" << endl;
     report( NOTICE, kModule ) << "*  due to random oversampling of the input file.   *" << endl;
     report( NOTICE, kModule ) << "*         Random Seed:  " << std::setw(7) << randSeed << "                    *" << endl;
     report( NOTICE, kModule ) << "****************************************************" << endl;
     report( NOTICE, kModule ) << endl;

      TH1::AddDirectory( kFALSE );
      gSystem->Load( "libTree" );

      // Open input file and tree
      ifstream fileexists( inFileName.c_str() );
      if (fileexists){
         fileexists.close();
         m_inFile = new TFile( inFileName.c_str() );
         if (!m_inFile || m_inFile->IsZombie()) {
            report( ERROR, kModule ) << "FSRootDataReaderBootstrap WARNING:  Cannot open file... " << inFileName << endl;
            m_inFile = NULL;
            m_inTree = NULL;
            return;
         }
         m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );

         if (!m_inTree) {
            report( ERROR, kModule ) << "FSRootDataReaderBootstrap WARNING:  Cannot open tree... " << inTreeName << endl;
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
         report( ERROR, kModule ) << "FSRootDataReaderBootstrap WARNING:  Cannot find file... " << inFileName << endl;
         m_inFile = NULL;
         m_inTree = NULL;
         return;
      }

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

      // Generate randomized event indices for bootstrap sampling
      unsigned int nEvents = numEvents();
      for( unsigned int i = 0; i < nEvents; ++i ){
         m_entryOrder.insert( (unsigned int)floor( m_randGenerator->Rndm() * nEvents ) );
      }

      m_nextEntry = m_entryOrder.begin(); // Initialize iterator
   }

void FSRootDataReaderBootstrap::resetSource(){
   m_eventCounter = 0;
   m_nextEntry = m_entryOrder.begin();  // Reset bootstrap iterator
}

Kinematics* FSRootDataReaderBootstrap::getEvent(){
   if( m_eventCounter < numEvents() ){

      // Modified: Use bootstrapped index instead of sequential order
      assert( m_nextEntry != m_entryOrder.end() );
      m_inTree->GetEntry( *m_nextEntry++ );

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

unsigned int FSRootDataReaderBootstrap::numEvents() const{
   if (!m_inTree) return 0;
   return static_cast< unsigned int >( m_inTree->GetEntries() );
}
