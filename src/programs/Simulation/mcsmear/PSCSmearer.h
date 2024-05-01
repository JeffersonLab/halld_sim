// Smearing class for coarse pair spectrometer counters (PSC)

#ifndef _PSCSMEARER_H_
#define _PSCSMEARER_H_

#include "Smearer.h"


class psc_config_t 
{
  public:
	psc_config_t(const std::shared_ptr<const JEvent>& event);

	double PSC_SIGMA;
	double PSC_PHOTONS_PERMEV;
	double PSC_THRESHOLD;

};


class PSCSmearer : public Smearer
{
  public:
	PSCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		psc_config = new psc_config_t(event);
	}
	~PSCSmearer() {
		delete psc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	psc_config_t  *psc_config;
};



#endif // _PSSMEARER_H_