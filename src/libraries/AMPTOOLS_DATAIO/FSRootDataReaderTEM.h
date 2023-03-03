#if !defined(FSROOTDATAREADERTEM)
#define FSROOTDATAREADERTEM

#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"

using namespace std;

class FSRootDataReaderTEM : public UserDataReader< FSRootDataReaderTEM >{

   public:

      FSRootDataReaderTEM() : UserDataReader< FSRootDataReaderTEM >() { }

      FSRootDataReaderTEM( const vector< string >& args );

      string name() const { return "FSRootDataReaderTEM"; }

      virtual Kinematics* getEvent();

      virtual void resetSource();

      virtual unsigned int numEvents() const;

      unsigned int eventCounter() const { return m_eventCounter; }


   private:

      	TFile* m_inFile;
      	TTree* m_inTree;
      	TTree* m_inFriendTree;
      	string inFileName,inTreeName;
	unsigned int m_eventCounter,m_numEvents;
      	unsigned int m_numParticles;
	double m_tMin, m_tMax, m_EMin, m_EMax, m_MMin, m_MMax;
	
      	double m_EnPB;
      	double m_PxPB;
      	double m_PyPB;
      	double m_PzPB;
      	double m_EnP[50];
      	double m_PxP[50];
      	double m_PyP[50];
      	double m_PzP[50];

      	double m_weight;

};

#endif
