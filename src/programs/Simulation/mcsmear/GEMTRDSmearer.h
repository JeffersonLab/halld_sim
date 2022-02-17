// Smearing class for GEMTRD

#ifndef _GEMTRDSMEARER_H_
#define _GEMTRDSMEARER_H_

#include "Smearer.h"

class gemtrd_config_t 
{
  public:
  gemtrd_config_t(JEventLoop *loop);
  
  // GEMTRD resolutions and threshold
  double GEMTRD_TSIGMA;
  double GEMTRD_ASIGMA;
  double GEMTRD_THRESHOLD;
  double GEMTRD_XYSIGMA;

  double GEMTRD_INTEGRAL_TO_AMPLITUDE;
};


class GEMTRDSmearer : public Smearer
{
  public:
	GEMTRDSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		gemtrd_config = new gemtrd_config_t(loop);
	}
	~GEMTRDSmearer() {
		delete gemtrd_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	gemtrd_config_t  *gemtrd_config;
};




#endif // _GEMTRDSMEARER_H_

