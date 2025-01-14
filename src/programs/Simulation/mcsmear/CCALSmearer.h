// Smearing class for compton calorimeter (CCAL)

#ifndef _CCALSMEARER_H_
#define _CCALSMEARER_H_

#include "Smearer.h"

#include <CCAL/DCCALGeometry.h>


class ccal_config_t 
{
  public:
	ccal_config_t(JEventLoop *loop);

	double CCAL_EN_SCALE;
	
	double CCAL_EN_P0;
	double CCAL_EN_P1;
	double CCAL_EN_P2;
	
	double CCAL_EN_GP0;
	double CCAL_EN_GP1;
	double CCAL_EN_GP2;
	
	vector<double> CCAL_PEDS;
	vector<double> CCAL_GAINS;
	
	double CCAL_THRESHOLD;
	double CCAL_INTEGRAL_PEAK;
	double CCAL_ADC_ASCALE;
	double CCAL_PED_RMS;
	
	// Time smearing factor
	double CCAL_TSIGMA;
	
	// Single block energy threshold (applied after smearing)
	double CCAL_BLOCK_THRESHOLD;
	
};


class CCALSmearer : public Smearer
{
  public:
	CCALSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		ccal_config = new ccal_config_t(loop);
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
