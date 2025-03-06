// Smearing class for compton calorimeter (CCAL)

#ifndef _CCALSMEARER_H_
#define _CCALSMEARER_H_

#include "Smearer.h"

#include <CCAL/DCCALGeometry.h>


class ccal_config_t 
{
  public:
	ccal_config_t(const std::shared_ptr<const JEvent>& event);

	double CCAL_EN_SCALE;
	
	double CCAL_EN_P0;
	double CCAL_EN_P1;
	double CCAL_EN_P2;	
	
	double CCAL_EN_GP0;
	double CCAL_EN_GP1;
	double CCAL_EN_GP2;

	
	// Time smearing factor
	double CCAL_TSIGMA;
	
	
	// Single block energy threshold (applied after smearing)
	double CCAL_BLOCK_THRESHOLD;
	
};


class CCALSmearer : public Smearer
{
  public:
	CCALSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		ccal_config = new ccal_config_t(event);
		ccalGeom = new DCCALGeometry();
	}
	~CCALSmearer() {
		delete ccal_config;
		delete ccalGeom;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	ccal_config_t  *ccal_config;
	DCCALGeometry  *ccalGeom;
};


#endif // _CCALSMEARER_H_
