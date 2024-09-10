// Smearing class for compton calorimeter (ECAL)

#ifndef _ECALSMEARER_H_
#define _ECALSMEARER_H_

#include "Smearer.h"

class ecal_config_t 
{
  public:
	ecal_config_t(JEventLoop *loop);

	double ECAL_EN_SCALE;
	
	double ECAL_EN_P0;
	double ECAL_EN_P1;
	double ECAL_EN_P2;	
	
	double ECAL_EN_GP0;
	double ECAL_EN_GP1;
	double ECAL_EN_GP2;

	
	// Time smearing factor
	double ECAL_TSIGMA;
	
	
	// Single block energy threshold (applied after smearing)
	double ECAL_BLOCK_THRESHOLD;
	
};


class ECALSmearer : public Smearer
{
  public:
        ECALSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
        ecal_config = new ecal_config_t(loop);
         }
	~ECALSmearer() {
		delete ecal_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	
  private:
  	ecal_config_t  *ecal_config;

};


#endif // _ECALSMEARER_H_
