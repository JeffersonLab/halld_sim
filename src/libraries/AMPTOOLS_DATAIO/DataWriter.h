#if !defined(DATAWRITER)
#define DATAWRITER

#include "IUAmpTools/Kinematics.h"

#include "TTree.h"
#include "TFile.h"

class DataWriter
{

public:

  DataWriter(){}
  virtual ~DataWriter(){}
  
  virtual void writeEvent( const Kinematics& kin ) = 0;
  virtual int eventCounter() const = 0;
  
private:

};

#endif
