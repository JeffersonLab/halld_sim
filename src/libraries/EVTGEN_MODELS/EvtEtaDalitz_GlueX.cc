// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>

#include "EvtEtaDalitz_GlueX.h"

EvtEtaDalitz_GlueX::~EvtEtaDalitz_GlueX() {}

std::string EvtEtaDalitz_GlueX::getName(){

  return "ETA_DALITZ_GLUEX";     

}


EvtDecayBase* EvtEtaDalitz_GlueX::clone(){

  return new EvtEtaDalitz_GlueX;

}

void EvtEtaDalitz_GlueX::init() 
{
	// check that there are 7 arguments
	checkNArg(7);
	checkNDaug(3);

	checkSpinParent(EvtSpinType::SCALAR);

	checkSpinDaughter(0,EvtSpinType::SCALAR);
	checkSpinDaughter(1,EvtSpinType::SCALAR);
	checkSpinDaughter(2,EvtSpinType::SCALAR);
}


void EvtEtaDalitz_GlueX::initProbMax()
{

	setProbMax(2.15);

}

void EvtEtaDalitz_GlueX::decay(EvtParticle *p)  
{

	p->initializePhaseSpace(getNDaug(), getDaugs());

	EvtVector4R mompip = p->getDaug(0)->getP4();
	EvtVector4R mompim = p->getDaug(1)->getP4();
	EvtVector4R mompi0 = p->getDaug(2)->getP4();
	double masspip = p->getDaug(0)->mass();
	double masspim = p->getDaug(1)->mass();
	double masspi0 = p->getDaug(2)->mass();
	double m_eta = p->mass();

	// Dalitz plot parameters
	double Q = m_eta-masspip-masspim-masspi0;
	double x = sqrt(3.)*( (mompip.get(0)-masspip) - (mompim.get(0)-masspim) )/Q; 
	double y = (mompi0.get(0)-masspi0)*(3.0/Q)-1.0;

	// pull out parameters from arguments
	double a = getArg(0);
	double b = getArg(1);
	double c = getArg(2);
	double d = getArg(3);
	double e = getArg(4);
	double f = getArg(5);
	double g = getArg(6);

	EvtComplex amp(sqrt(1.0 + a*y + b*y*y + c*x + d*x*x + e*x*y + f*y*y*y + g*x*x*y),0.0);

	vertex(amp);

	return;

}


