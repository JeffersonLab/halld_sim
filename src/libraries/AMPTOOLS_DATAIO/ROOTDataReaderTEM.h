#if !defined(ROOTDATAREADERTEM)
#define ROOTDATAREADERTEM

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReaderTEM : public UserDataReader< ROOTDataReaderTEM >
{
	
public:
  
  /**
   * Default constructor for ROOTDataReaderTEM
   */
  ROOTDataReaderTEM() : UserDataReader< ROOTDataReaderTEM >(), m_inFile( NULL ) { }
  
  ~ROOTDataReaderTEM();
  
  /**
   * Constructor for ROOTDataReaderTEM
   * \param[in] args vector of string arguments
   */
  ROOTDataReaderTEM( const vector< string >& args );
  
  string name() const { return "ROOTDataReaderTEM"; }
 
  virtual vector<TLorentzVector> particleList();
  virtual bool checkEvent();
  virtual Kinematics* getEvent();
  virtual void resetSource();

  /**
   * This function returns a true if the file was open
   * with weight-reading enabled and had this tree branch,
   * false, if these criteria are not met.
   */
  virtual bool hasWeight(){ return m_useWeight; };
  virtual unsigned int numEvents() const;
  
private:
	
  TFile* m_inFile;
  TTree* m_inTree;
  unsigned int m_eventCounter,m_numEvents;
  bool m_useWeight, m_RangeSpecified;
  double m_tMin,m_tMax, m_EMin,m_EMax, m_MMin,m_MMax;
  
  int m_nPart;
  float m_e[Kinematics::kMaxParticles];
  float m_px[Kinematics::kMaxParticles];
  float m_py[Kinematics::kMaxParticles];
  float m_pz[Kinematics::kMaxParticles];
  float m_eBeam;
  float m_pxBeam;
  float m_pyBeam;
  float m_pzBeam;
  float m_weight;
};

#endif
