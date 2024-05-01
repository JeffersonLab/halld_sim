// Smearing class for central drift chamber (CDC)

#ifndef _CDCSMEARER_H_
#define _CDCSMEARER_H_

#include "Smearer.h"


class cdc_config_t 
{
  public:
	cdc_config_t(const std::shared_ptr<const JEvent>& event);

	double CDC_TDRIFT_SIGMA;
	double CDC_TIME_WINDOW;
	double CDC_PEDESTAL_SIGMA;   // deprecated
	double CDC_THRESHOLD_FACTOR; // number of pedestal sigmas for determining sparsification threshold - deprecated
        double CDC_INTEGRAL_TO_AMPLITUDE;
        double CDC_ASCALE;
	double CDC_DIFFUSION_PAR1,CDC_DIFFUSION_PAR2,CDC_DIFFUSION_PAR3;

        vector<double> CDC_GAIN_DOCA_PARS;  // params to model gas deterioration spring 2018
        vector<double> CDC_GAIN_DOCA_EXT;  // params to model gas deterioration spring 2018

	vector< vector<double> > wire_efficiencies;
	vector< vector<double> > wire_thresholds;

	void CalcNstraws(const std::shared_ptr<const JEvent>& event, int32_t runnumber, vector<unsigned int> &Nstraws);
	double GetEfficiencyCorrectionFactor(int ring, int straw) {
		return wire_efficiencies.at(ring-1).at(straw-1);
	}
	double GetWireThreshold(int ring, int straw) {
		return wire_thresholds.at(ring-1).at(straw-1);
	}
};


class CDCSmearer : public Smearer
{
  public:
	CDCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		cdc_config = new cdc_config_t(event);
	}
	~CDCSmearer() {
		delete cdc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	cdc_config_t  *cdc_config;
};



#endif // _CDCSMEARER_H_
