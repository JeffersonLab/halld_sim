// Smearing class for forward calorimeter (CGEM)

#ifndef _CGEMSMEARER_H_
#define _CGEMSMEARER_H_

#include "Smearer.h"


class cgem_config_t 
{
  public:
	cgem_config_t(JEventLoop *loop);
	
};


class CGEMSmearer : public Smearer
{
 public:
 CGEMSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
    cgem_config = new cgem_config_t(loop);
  }
  ~CGEMSmearer() {
    delete cgem_config;
  }
  
  void SmearEvent(hddm_s::HDDM *record);
  
 private:
  cgem_config_t *cgem_config;

};


#endif // _CGEMSMEARER_H_
