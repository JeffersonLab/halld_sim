// Smearing class for tagger microscope (TAGM)

#ifndef _TAGMSMEARER_H_
#define _TAGMSMEARER_H_

#include "Smearer.h"


class tagm_config_t 
{
  public:
	tagm_config_t(const std::shared_ptr<const JEvent>& event);

	double TAGM_TSIGMA;
	double TAGM_FADC_TSIGMA;
	double TAGM_NPIX_PER_GEV;

    std::map<int, int> fiber_quality;
    std::map<int, std::vector<double> > energy_range_GeV;
    double endpoint_energy_GeV;
    double endpoint_calib_GeV;
};


class TAGMSmearer : public Smearer
{
  public:
	TAGMSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		tagm_config = new tagm_config_t(event);
	}
	~TAGMSmearer() {
		delete tagm_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

    static double get_tagm_energy(int column, int low_mid_high=-1);

  private:
  	tagm_config_t  *tagm_config;
};



#endif // _TAGMSMEARER_H_
