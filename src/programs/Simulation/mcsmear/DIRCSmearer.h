// Smearing class for DIRC

#ifndef _DIRCSMEARER_H_
#define _DIRCSMEARER_H_

#include "Smearer.h"

class dirc_config_t
{
  public:
        dirc_config_t(JEventLoop *loop);

        double DIRC_TSIGMA;
	double DIRC_EFFIC_SCALE;

	int DIRC_MAX_CHANNELS;
	vector< vector <int> > dChannelStatus;
	vector< vector <float> > dChannelEffic;
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
	enum dirc_status_state {GOOD, BAD, NOISY};
};


#endif // _DIRCSMEARER_H_
