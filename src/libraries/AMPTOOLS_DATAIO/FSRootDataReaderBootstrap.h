#if !defined(FSROOTDATAREADERBOOTSTRAP)
#define FSROOTDATAREADERBOOTSTRAP

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

class FSRootDataReaderBootstrap : public UserDataReader< FSRootDataReaderBootstrap >{

   public:

      FSRootDataReaderBootstrap() : UserDataReader< FSRootDataReaderBootstrap >() { }

      FSRootDataReaderBootstrap( const vector< string >& args );

      string name() const { return "FSRootDataReaderBootstrap"; }

      virtual Kinematics* getEvent();

      virtual void resetSource();

      virtual unsigned int numEvents() const;

      unsigned int eventCounter() const { return m_eventCounter; }

   private:

      TFile* m_inFile;
      TTree* m_inTree;
      TTree* m_inFriendTree;
      unsigned int m_eventCounter;
      unsigned int m_numParticles;

      double m_EnPB;
      double m_PxPB;
      double m_PyPB;
      double m_PzPB;
      double m_EnP[50];
      double m_PxP[50];
      double m_PyP[50];
      double m_PzP[50];

      double m_weight;

      // Added for bootstrap functionality
      TRandom3* m_randGenerator;  // Random number generator for bootstrapping
      std::multiset<unsigned int> m_entryOrder;  // Stores bootstrap sampled event indices
      mutable std::multiset<unsigned int>::const_iterator m_nextEntry;  // Iterator for sampling

};

#endif
