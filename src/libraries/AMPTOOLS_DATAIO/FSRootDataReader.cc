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

FSRootDataReader::FSRootDataReader( const vector< string >& args ) :
   UserDataReader< FSRootDataReader >(args),
   m_eventCounter( 0 ){

      assert((args.size() >= 3 && args.size() <= 4) || (args.size()>=6 && args.size()<=7));
      string inFileName(args[0]);
      string inTreeName(args[1]);
      m_numParticles = atoi(args[2].c_str());
      assert (m_numParticles < 50);
      TString fourMomentumPrefix = "";
      if (args.size() == 4) fourMomentumPrefix = args[3];
      
      TString friendFileName = "";
      TString friendBranchName = "weight";
      TString friendTreeName = "";
      if (args.size() >= 6) {
        friendFileName = args[3];
        friendTreeName = args[4];
        friendBranchName = args[5];
      }
      if (args.size() == 7) fourMomentumPrefix = args[6];

      TH1::AddDirectory( kFALSE );
      gSystem->Load( "libTree" );

      ifstream fileexists( inFileName.c_str() );
      if (fileexists){
         m_inFile = new TFile( inFileName.c_str() );
         m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
         if(args.size()>=5)
            m_inTree->AddFriend(friendTreeName, friendFileName);
      }
      else{
         cout << "FSRootDataReader WARNING:  Cannot find file... " << inFileName << endl;
         m_inFile = NULL;
         m_inTree = NULL;
      }

      if(args.size()==3)
        cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << endl;
      if(args.size()==4)
        cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << " " << args[3] << endl;
      if(args.size()==5)
        cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << " " << args[3] << " " << args[4] << endl;
      if(args.size()==6)
        cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << " " << args[3] << " " << args[4] << " " << args[5] << endl;
      if(args.size()==7)
        cout << "Opening Tree " << args[0] << " " << args[1] << " " << args[2] << " " << args[3] << " " << args[4] << " " << args[5] << " " << args[6] << endl;
      
      if (m_inTree){

         TString sEnPB = fourMomentumPrefix+"EnPB";
         TString sPxPB = fourMomentumPrefix+"PxPB";
         TString sPyPB = fourMomentumPrefix+"PyPB";
         TString sPzPB = fourMomentumPrefix+"PzPB";
         m_inTree->SetBranchAddress( sEnPB, &m_EnPB );
         m_inTree->SetBranchAddress( sPxPB, &m_PxPB );
         m_inTree->SetBranchAddress( sPyPB, &m_PyPB );
         m_inTree->SetBranchAddress( sPzPB, &m_PzPB );

         // particle name labels, matching ROOT branch convention
         std::vector<TString> particleLabels;

         if (m_numParticles == 5) {
            particleLabels = { 
                "1",   // Lambda
                "2",   // K+
                "3",   // Pi0
                "1a",  // Proton from Lambda decay
                "1b"   // Pi- from Lambda decay
            };
         } 
         else if (m_numParticles == 3) {
               particleLabels = { 
                  "2",    // K+
                  "1a",  // Proton from Lambda decay
                  "1b"  // Pi- from Lambda decay
               };
         }

         for (unsigned int i = 0; i < m_numParticles; i++) {
            TString label = particleLabels[i];
            TString sEnPi = fourMomentumPrefix + "EnP" + label;
            TString sPxPi = fourMomentumPrefix + "PxP" + label;
            TString sPyPi = fourMomentumPrefix + "PyP" + label;
            TString sPzPi = fourMomentumPrefix + "PzP" + label;
        
            m_inTree->SetBranchAddress(sEnPi, &m_EnP[i]);
            m_inTree->SetBranchAddress(sPxPi, &m_PxP[i]);
            m_inTree->SetBranchAddress(sPyPi, &m_PyP[i]);
            m_inTree->SetBranchAddress(sPzPi, &m_PzP[i]);
        }
        
        // handle event weight
        if (args.size() >= 6)
            m_inTree->SetBranchAddress(friendBranchName, &m_weight);
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
//      m_weight = 1.0;
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
