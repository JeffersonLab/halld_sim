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

      // these are default initialized so that a four-momentum which never gets
      // bound to a branch reads back as zero rather than as uninitialized memory
      TFile* m_inFile = nullptr;
      TTree* m_inTree = nullptr;
      TTree* m_inFriendTree = nullptr;
      unsigned int m_eventCounter = 0;
      unsigned int m_numParticles = 0;

      double m_EnPB = 0.;
      double m_PxPB = 0.;
      double m_PyPB = 0.;
      double m_PzPB = 0.;
      double m_EnP[50] = {};
      double m_PxP[50] = {};
      double m_PyP[50] = {};
      double m_PzP[50] = {};

      double m_weight = 1.;

      // Added for bootstrap functionality
      TRandom3* m_randGenerator = nullptr;  // Random number generator for bootstrapping
      std::multiset<unsigned int> m_entryOrder;  // Stores bootstrap sampled event indices
      mutable std::multiset<unsigned int>::const_iterator m_nextEntry;  // Iterator for sampling

      static const char* kModule;
};

#endif
