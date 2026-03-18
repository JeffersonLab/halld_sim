// Smearing class for forward multiple-wire proportional chamber (FMWPC)

#ifndef _FMWPCSMEARER_H_
#define _FMWPCSMEARER_H_

#include "Smearer.h"

class fmwpc_config_t 
{
  public:
  fmwpc_config_t(const std::shared_ptr<const JEvent>& event);
  
  // FMWPC resolutions and threshold
  double FMWPC_TSIGMA;
  double FMWPC_ASIGMA;
  double FMWPC_THRESHOLD;
};


class FMWPCSmearer : public Smearer
{
  public:
	FMWPCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		fmwpc_config = new fmwpc_config_t(event);
	}
	~FMWPCSmearer() {
		delete fmwpc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	fmwpc_config_t  *fmwpc_config;
};




#endif // _FMWPCSMEARER_H_

