#if !defined(ROOTDATAREADERHIST)
#define ROOTDATAREADERHIST

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <string>

using namespace std;

class ROOTDataReaderHist : public UserDataReader< ROOTDataReaderHist >
{
	
public:
  
  /**
   * Default constructor for ROOTDataReaderHist
   */
  ROOTDataReaderHist() : UserDataReader< ROOTDataReaderHist >(), m_inFile( NULL ) { }
  
  ~ROOTDataReaderHist();
  
  /**
   * Constructor for ROOTDataReaderHist
   * \param[in] args vector of string arguments
   */
  ROOTDataReaderHist( const vector< string >& args );
  
  string name() const { return "ROOTDataReaderHist"; }
  
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
  unsigned int m_eventCounter;
  bool m_useWeight;
  
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

  TFile* m_HistFile;
  TH1F* m_Hist;
};

#endif
