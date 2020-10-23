
#ifndef EVTETADALITZ_3PI0_GLUEX_HH
#define EVTETADALITZ_3PI0_GLUEX_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtEtaDalitz_3Pi0_GlueX:public  EvtDecayAmp  {

public:

  EvtEtaDalitz_3Pi0_GlueX() {}
  virtual ~EvtEtaDalitz_3Pi0_GlueX();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 
  
private:
	
	int warning_counter;

};

#endif
