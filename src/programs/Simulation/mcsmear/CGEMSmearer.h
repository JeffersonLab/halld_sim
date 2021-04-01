// Smearing class for forward calorimeter (CGEM)

#ifndef _CGEMSMEARER_H_
#define _CGEMSMEARER_H_

#include "Smearer.h"


class cgem_config_t 
{
  public:
	cgem_config_t(JEventLoop *loop);
		

	double m_CGEM_WORK_FUNCTION;
	double m_CGEM_FANO_FACTOR; 
	double m_CGEM_ENERGY_THRES;
	double m_CGEM_SPATIAL_RESOLUTION;
	double m_CGEM_Z_RESOLUTION;
	double m_CGEM_TIMING_RESOLUTION;

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
