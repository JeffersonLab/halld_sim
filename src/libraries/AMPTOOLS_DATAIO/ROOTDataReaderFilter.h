#if !defined(ROOTDATAREADERFILTER)
#define ROOTDATAREADERFILTER

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReaderFilter : public UserDataReader< ROOTDataReaderFilter >
{
	
public:
  
  /**
   * Default constructor for ROOTDataReaderFilter
   */
  ROOTDataReaderFilter() : UserDataReader< ROOTDataReaderFilter >(), m_inFile( NULL ) { }
  
  ~ROOTDataReaderFilter();
  
  /**
   * Constructor for ROOTDataReaderFilter
   * \param[in] args vector of string arguments
   */
  ROOTDataReaderFilter( const vector< string >& args );
  
  string name() const { return "ROOTDataReaderFilter"; }
  
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
  double m_tMin,m_tMax;

  int nargs;
  static const int maxselects=25; // maximum number of selections you wish to apply
  int nselects;
  map<string,float> mapVars={};
  float s_min[maxselects];
  float s_max[maxselects];
  bool b_isSelection[maxselects];
  string s_var[maxselects];
  bool selection;
  
  int m_nPart;
  float m_weightIntegral;
  float m_e[Kinematics::kMaxParticles];
  float m_px[Kinematics::kMaxParticles];
  float m_py[Kinematics::kMaxParticles];
  float m_pz[Kinematics::kMaxParticles];
  float m_eBeam;
  float m_pxBeam;
  float m_pyBeam;
  float m_pzBeam;
  float m_weight=1;
};

#endif