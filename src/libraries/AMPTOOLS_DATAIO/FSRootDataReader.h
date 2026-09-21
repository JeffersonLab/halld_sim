#if !defined(FSROOTDATAREADER)
#define FSROOTDATAREADER

#include <string>
#include <set>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom2.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"

using namespace std;

class FSRootDataReader : public UserDataReader< FSRootDataReader >{

   public:

      FSRootDataReader() : UserDataReader< FSRootDataReader >() { }

      FSRootDataReader( const vector< string >& args );

      string name() const { return "FSRootDataReader"; }

      virtual Kinematics* getEvent();

      virtual void resetSource();

      virtual void resample( unsigned int seed );

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

      static const char* kModule;
      TRandom2* m_randGenerator = NULL;
      multiset< unsigned int > m_entryOrder;
      mutable multiset< unsigned int >::const_iterator m_nextEntry;
};

#endif
