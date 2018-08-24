// Smearing class for DIRC

#ifndef _DIRCSMEARER_H_
#define _DIRCSMEARER_H_

#include "Smearer.h"

#include <DIRC/DDIRCGeometry.h>

class dirc_config_t
{
  public:
        dirc_config_t(JEventLoop *loop, DDIRCGeometry *dircGeom);

        double DIRC_TSIGMA;
};


class DIRCSmearer : public Smearer
{
  public:
	DIRCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		dircGeom = new DDIRCGeometry(loop->GetJEvent().GetRunNumber());
                dirc_config = new dirc_config_t(loop, dircGeom);
	}
	~DIRCSmearer() {
		delete dircGeom;
		delete dirc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

 private:
        dirc_config_t  *dirc_config;
        DDIRCGeometry *dircGeom;
};


#endif // _DIRCSMEARER_H_
