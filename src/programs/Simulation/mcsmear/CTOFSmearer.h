// Smearing class for CPP scintillator paddles

#ifndef _CTOFSMEARER_H_
#define _CTOFSMEARER_H_

#include "Smearer.h"

#include <TMath.h>

class ctof_config_t 
{
  public:
  ctof_config_t(const std::shared_ptr<const JEvent>& event);

  double TSIGMA;
  double PHOTONS_PERMEV;
  double BAR_THRESHOLD;
  double ATTENUATION_LENGTH;
  double BAR_LENGTH;
};


class CTOFSmearer : public Smearer
{
  public:
 CTOFSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
    ctof_config = new ctof_config_t(event);
  }
  ~CTOFSmearer() {
    delete ctof_config;
  }
  
  void SmearEvent(hddm_s::HDDM *record);
  
 private:
  ctof_config_t  *ctof_config;
  
};


#endif // _CTOFSMEARER_H_
