// Smearing class for DIRC

#ifndef _DIRCSMEARER_H_
#define _DIRCSMEARER_H_

#include "Smearer.h"

class dirc_config_t
{
  public:
        dirc_config_t(JEventLoop *loop);

        double DIRC_TSIGMA;
};


class DIRCSmearer : public Smearer
{
  public:
	DIRCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
                dirc_config = new dirc_config_t(loop);
	}
	~DIRCSmearer() {
		delete dirc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

 private:
        dirc_config_t  *dirc_config;
};


#endif // _DIRCSMEARER_H_
