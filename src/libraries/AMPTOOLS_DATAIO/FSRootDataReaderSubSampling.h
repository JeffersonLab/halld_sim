#if !defined(FSROOTDATAREADERSUBSAMPLING)
#define FSROOTDATAREADERSUBSAMPLING

#include <string>
#include <vector>
#include <set>  // Added for multiset
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"  // Added for random number generator
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"

using namespace std;

class FSRootDataReaderSubSampling : public UserDataReader< FSRootDataReaderSubSampling >{

   public:

      FSRootDataReaderSubSampling() : UserDataReader< FSRootDataReaderSubSampling >() { }

      FSRootDataReaderSubSampling( const vector< string >& args );

      string name() const { return "FSRootDataReaderSubSampling"; }

      Kinematics* getEvent() override;

      void resetSource() override;

      unsigned int numEvents() const override;

      unsigned int eventCounter() const { return m_eventCounter; }

   private:

      TFile* m_inFile;
      TTree* m_inTree;
      TTree* m_inFriendTree;
      unsigned int m_eventCounter;
      unsigned int m_numParticles;
      double m_subsampleFraction;

      double m_EnPB;
      double m_PxPB;
      double m_PyPB;
      double m_PzPB;
      double m_EnP[50];
      double m_PxP[50];
      double m_PyP[50];
      double m_PzP[50];

      double m_weight;

      // Added for SubSampling functionality
      TRandom3* m_randGenerator;  // Random number generator for SubSamplingping
      // std::multiset<unsigned int> m_entryOrder;  // Stores SubSampling sampled event indices
      // mutable std::multiset<unsigned int>::const_iterator m_nextEntry;  // Iterator for sampling
      std::vector<unsigned int> m_entryOrder;  // Stores sampled event indices
      mutable std::vector<unsigned int>::const_iterator m_nextEntry;  // Iterator for sampled entries

};

#endif
