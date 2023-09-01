// Smearing class for inner TOF detector

#ifndef _ITOFSMEARER_H_
#define _ITOFSMEARER_H_

#include "Smearer.h"

#include <TMath.h>

class itof_config_t 
{
  public:
  itof_config_t(JEventLoop *loop);

  double TSIGMA;
  double PHOTONS_PERMEV;
};


class ITOFSmearer : public Smearer
{
  public:
 ITOFSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
    itof_config = new itof_config_t(loop);
  }
  ~ITOFSmearer() {
    delete itof_config;
  }
  
  void SmearEvent(hddm_s::HDDM *record);
  
 private:
  itof_config_t  *itof_config;
  
};


#endif // _ITOFSMEARER_H_
