// Smearing class for tagger hodoscope (TAGH)

#ifndef _TAGHSMEARER_H_
#define _TAGHSMEARER_H_

#include "Smearer.h"


class tagh_config_t 
{
  public:
	tagh_config_t(JEventLoop *loop);

	double TAGH_TSIGMA;
	double TAGH_FADC_TSIGMA;
	double TAGH_NPE_PER_GEV;

    std::map<int, int> counter_quality;
    std::map<int, std::vector<double> > energy_range_GeV;
    double endpoint_energy_GeV;
    double endpoint_calib_GeV;
};


class TAGHSmearer : public Smearer
{
  public:
	TAGHSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		tagh_config = new tagh_config_t(loop);
	}
	~TAGHSmearer() {
		delete tagh_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

    static double get_tagh_energy(int counter, int low_mid_high=-1);

  private:
  	tagh_config_t  *tagh_config;
};




#endif // _TAGHSMEARER_H_
