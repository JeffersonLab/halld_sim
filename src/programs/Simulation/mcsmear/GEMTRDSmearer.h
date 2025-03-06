// Smearing class for GEMTRD

#ifndef _GEMTRDSMEARER_H_
#define _GEMTRDSMEARER_H_

#include "Smearer.h"

class gemtrd_config_t 
{
  public:
  gemtrd_config_t(const std::shared_ptr<const JEvent>& event);
  
  // GEMTRD resolutions and threshold
  double GEMTRD_TSIGMA;
  double GEMTRD_ASIGMA;
  double GEMTRD_THRESHOLD;
  double GEMTRD_XYSIGMA;

  double GEMTRD_INTEGRAL_TO_AMPLITUDE;
};


class GEMTRDSmearer : public Smearer
{
public:
  GEMTRDSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t 
*in_config) : Smearer(event, in_config) {
    gemtrd_config = new gemtrd_config_t(event);
  }

  ~GEMTRDSmearer() {
    delete gemtrd_config;
  }
  
  void SmearEvent(hddm_s::HDDM *record);
  double GetStripCharge(double q,int strip) const;
  
private:
  class strip_hit_t{
  public:
    strip_hit_t(int plane,int strip,double q,double t):plane(plane),strip(strip),q(q),t(t){}
    int plane;
    int strip;
    double q;
    double t;
  };

  gemtrd_config_t  *gemtrd_config;
};




#endif // _GEMTRDSMEARER_H_

