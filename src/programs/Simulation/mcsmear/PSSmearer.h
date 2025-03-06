// Smearing class for fine pair spectrometer counters (PS)

#ifndef _PSSMEARER_H_
#define _PSSMEARER_H_

#include "Smearer.h"


class ps_config_t 
{
  public:
	ps_config_t(const std::shared_ptr<const JEvent>& event);

	double PS_SIGMA;
	double PS_NPIX_PER_GEV;
	double PS_THRESHOLD;

};


class PSSmearer : public Smearer
{
  public:
	PSSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		ps_config = new ps_config_t(event);
	}
	~PSSmearer() {
		delete ps_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	ps_config_t  *ps_config;
};


#endif // _PSSMEARER_H_