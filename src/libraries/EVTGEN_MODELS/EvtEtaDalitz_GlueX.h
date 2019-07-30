
#ifndef EVTETADALITZ_GLUEX_HH
#define EVTETADALITZ_GLUEX_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtEtaDalitz_GlueX:public  EvtDecayAmp  {

public:

  EvtEtaDalitz_GlueX() {}
  virtual ~EvtEtaDalitz_GlueX();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

};

#endif
