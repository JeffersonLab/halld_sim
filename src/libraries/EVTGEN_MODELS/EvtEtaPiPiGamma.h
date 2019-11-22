
#ifndef EVTETAPIPIGAMMA_HH
#define EVTETAPIPIGAMMA_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtEtaPiPiGamma:public  EvtDecayAmp  {

public:

  EvtEtaPiPiGamma() {}
  virtual ~EvtEtaPiPiGamma();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

};

#endif
