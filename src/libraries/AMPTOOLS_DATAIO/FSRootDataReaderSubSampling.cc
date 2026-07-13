#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <set>  // Added for multiset
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/FSRootDataReaderSubSampling.h"
#include "TSystem.h"
using namespace std;


FSRootDataReaderSubSampling::FSRootDataReaderSubSampling( const vector< string >& args ) :
   UserDataReader< FSRootDataReaderSubSampling >(args),
   m_eventCounter( 0 ){

      // args[0] = subsampleFraction
      assert((args.size() >= 4 && args.size() <= 6) || (args.size() >= 7 && args.size() <= 9));

      m_subsampleFraction = std::stod(args[0]);
      assert(m_subsampleFraction > 0.0 && m_subsampleFraction < 1.0);

      string inFileName(args[1]);
      string inTreeName(args[2]);
      m_numParticles = atoi(args[3].c_str());
      assert(m_numParticles < 50);

      TString fourMomentumPrefix = "";
      TString friendFileName = "";
      TString friendTreeName = "";
      TString friendBranchName = "weight";
      int randSeed = 0;

      auto isInteger = [](const string& s) {
         return !s.empty() && s.find_first_not_of("-0123456789") == string::npos;
      };

      // ---- shifted logic ----
      if (args.size() == 5) {
         if (isInteger(args[4])) {
            randSeed = atoi(args[4].c_str());
         } else {
            fourMomentumPrefix = args[4];
         }
      } 
      else if (args.size() == 6) {
         fourMomentumPrefix = args[4];
         randSeed = atoi(args[5].c_str());
      } 
      else if (args.size() >= 7) {
         friendFileName = args[4];
         friendTreeName = args[5];
         friendBranchName = args[6];

         if (args.size() == 8) {
            fourMomentumPrefix = args[7];
         }
         else if (args.size() == 9) {
            fourMomentumPrefix = args[7];
            randSeed = atoi(args[8].c_str());
         }
      }

      // Initialize random number generator with provided seed
      m_randGenerator = new TRandom3(randSeed);
      std::mt19937 rng(randSeed);

      // Compact printing of parsed arguments
      std::cout << "FSRootDataReaderSubSampling initialized with:\n"
          << "  inFileName:        " << inFileName << "\n"
          << "  inTreeName:        " << inTreeName << "\n"
          << "  numParticles:      " << m_numParticles << "\n"
          << "  friendFileName:    " << (friendFileName.Length() == 0 ? "N/A" : friendFileName.Data()) << "\n"
          << "  friendTreeName:    " << (friendTreeName.Length() == 0 ? "N/A" : friendTreeName.Data()) << "\n"
          << "  friendBranchName:  " << (friendBranchName.Length() == 0 ? "N/A" : friendBranchName.Data()) << "\n"
          << "  fourMomentumPrefix:" << (fourMomentumPrefix.Length() == 0 ? "N/A" : fourMomentumPrefix.Data()) << "\n"
          << "  randSeed:          " << randSeed << "\n"
          << "  sampling Fraction: " << m_subsampleFraction << "\n"
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
         if(args.size() >= 6)
            m_inTree->AddFriend(friendTreeName, friendFileName);
      }
      else{
         cout << "FSRootDataReaderSubSampling WARNING:  Cannot find file... " << inFileName << endl;
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
        if (args.size() >= 7)
            m_inTree->SetBranchAddress(friendBranchName, &m_weight);
        else
            m_weight = 1.0;
      }

      // Generate a deterministic random subsample without replacement
      if (!(m_subsampleFraction > 0.0 && m_subsampleFraction <= 1.0)) {
         std::cerr << "FSRootDataReaderSubSampling ERROR: subsample fraction must be in (0, 1], got "
                  << m_subsampleFraction << std::endl;
         std::exit(EXIT_FAILURE);
      }

      unsigned int nEvents = numEvents();
      m_entryOrder.clear(); m_entryOrder.reserve(nEvents);

      for (unsigned int i = 0; i < nEvents; ++i) {
         m_entryOrder.push_back(i);
      }

      // only sample fraction of total event after shuffled
      unsigned int nKeep = static_cast<unsigned int>(std::floor(m_subsampleFraction * nEvents));
      std::shuffle(m_entryOrder.begin(), m_entryOrder.end(), rng);
      m_entryOrder.resize(nKeep);

      // Sort indices for sequential TTree access to improve cache locality and I/O efficiency
      std::sort(m_entryOrder.begin(), m_entryOrder.end());
      m_nextEntry = m_entryOrder.begin();
   }

void FSRootDataReaderSubSampling::resetSource(){
   m_eventCounter = 0;
   m_nextEntry = m_entryOrder.begin();  // Reset SubSampling iterator
}

Kinematics* FSRootDataReaderSubSampling::getEvent() {
   if (m_nextEntry == m_entryOrder.end() || m_eventCounter >= m_entryOrder.size()) {
      return NULL;
   }

   m_inTree->GetEntry(*m_nextEntry++);
   ++m_eventCounter;

   vector<TLorentzVector> particleList;
   particleList.push_back(TLorentzVector(m_PxPB, m_PyPB, m_PzPB, m_EnPB));

   for (unsigned int i = 0; i < m_numParticles; i++) {
      particleList.push_back(TLorentzVector(m_PxP[i], m_PyP[i], m_PzP[i], m_EnP[i]));
   }

   return new Kinematics(particleList, m_weight);
}

unsigned int FSRootDataReaderSubSampling::numEvents() const {
    if (!m_inTree) return 0;

    // Before subsampling is initialized, fall back to full tree size
    if (m_entryOrder.empty()) {
        return static_cast<unsigned int>(m_inTree->GetEntries());
    }

    return static_cast<unsigned int>(m_entryOrder.size());
}
