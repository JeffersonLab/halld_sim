// Smearing class for DIRC

#ifndef _DIRCSMEARER_H_
#define _DIRCSMEARER_H_

#include "Smearer.h"

class dirc_config_t
{
  public:
        dirc_config_t(const std::shared_ptr<const JEvent>& event);

        double DIRC_TSIGMA;

	int DIRC_MAX_CHANNELS;
	vector< vector <int> > dChannelStatus;
};


class DIRCSmearer : public Smearer
{
  public:
	DIRCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
                dirc_config = new dirc_config_t(event);
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
