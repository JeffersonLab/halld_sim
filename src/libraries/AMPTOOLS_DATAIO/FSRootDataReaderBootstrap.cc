#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>  // Added for multiset
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/FSRootDataReaderBootstrap.h"
#include "TSystem.h"

using namespace std;

#include <iostream> // Include for printing

FSRootDataReaderBootstrap::FSRootDataReaderBootstrap( const vector< string >& args ) :
   UserDataReader< FSRootDataReaderBootstrap >(args),
   m_eventCounter( 0 ){

      assert((args.size() >= 3 && args.size() <= 5) || (args.size() >= 6 && args.size() <= 8));
      string inFileName(args[0]);
      string inTreeName(args[1]);
      m_numParticles = atoi(args[2].c_str());
      assert (m_numParticles < 50);
      
      TString fourMomentumPrefix = "";
      TString friendFileName = "";
      TString friendTreeName = "";
      TString friendBranchName = "weight";
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
          friendBranchName = args[5];

          if (args.size() == 7) {
              fourMomentumPrefix = args[6];
          }
          else if (args.size() == 8) {
              fourMomentumPrefix = args[6];
              randSeed = atoi(args[7].c_str());
          }
      }

      // Initialize random number generator with provided seed
      m_randGenerator = new TRandom3(randSeed);

      // Compact printing of parsed arguments
      std::cout << "FSRootDataReaderBootstrap initialized with:\n"
          << "  inFileName:        " << inFileName << "\n"
          << "  inTreeName:        " << inTreeName << "\n"
          << "  numParticles:      " << m_numParticles << "\n"
          << "  friendFileName:    " << (friendFileName.Length() == 0 ? "N/A" : friendFileName.Data()) << "\n"
          << "  friendTreeName:    " << (friendTreeName.Length() == 0 ? "N/A" : friendTreeName.Data()) << "\n"
          << "  friendBranchName:  " << (friendBranchName.Length() == 0 ? "N/A" : friendBranchName.Data()) << "\n"
          << "  fourMomentumPrefix:" << (fourMomentumPrefix.Length() == 0 ? "N/A" : fourMomentumPrefix.Data()) << "\n"
          << "  randSeed:          " << randSeed << "\n"
          << std::endl;


      cout << "******************** WARNING ***********************" << endl;
      cout << "*  You are using the boostrap data reader, which   *" << endl;
      cout << "*  should only be used for evaluating errors.      *" << endl;
      cout << "*  The results with different seeds will be random *" << endl;
      cout << "*  due to random oversampling of the input file.   *" << endl;
      cout << "****************************************************" << endl;
      cout << endl;
      cout << "   Random Seed:  " << randSeed << endl << endl;


      TH1::AddDirectory( kFALSE );
      gSystem->Load( "libTree" );

      ifstream fileexists( inFileName.c_str() );
      if (fileexists){
         m_inFile = new TFile( inFileName.c_str() );
         m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
         if(args.size() >= 5)
            m_inTree->AddFriend(friendTreeName, friendFileName);
      }
      else{
         cout << "FSRootDataReaderBootstrap WARNING:  Cannot find file... " << inFileName << endl;
         m_inFile = NULL;
         m_inTree = NULL;
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
            if(args.size() >= 6)
              m_inTree->SetBranchAddress( friendBranchName, &m_weight );
            else
              m_weight = 1.0;
         }
      }

      // Added: Generate randomized event indices for bootstrap sampling
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
