// Smearing class for TRD

#ifndef _TRDSMEARER_H_
#define _TRDSMEARER_H_

#include "Smearer.h"

class trd_config_t 
{
  public:
  trd_config_t(JEventLoop *loop);
  
  // TRD resolutions and threshold
  double TRD_TSIGMA;
  double TRD_ASIGMA;
  double TRD_THRESHOLD;
  double TRD_XYSIGMA;

  double TRD_INTEGRAL_TO_AMPLITUDE;
};


class TRDSmearer : public Smearer
{
  public:
	TRDSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		trd_config = new trd_config_t(loop);
	}
	~TRDSmearer() {
		delete trd_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	trd_config_t  *trd_config;
};




#endif // _TRDSMEARER_H_

