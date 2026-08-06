#if !defined(FSROOTDATAREADER)
#define FSROOTDATAREADER

#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
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

      static const char* kModule;
};

#endif
