// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>

#include "EvtEtaPiPiGamma.h"

EvtEtaPiPiGamma::~EvtEtaPiPiGamma() {}

std::string EvtEtaPiPiGamma::getName(){

  return "ETA_PIPIGAMMA";     

}


EvtDecayBase* EvtEtaPiPiGamma::clone(){

  return new EvtEtaPiPiGamma;

}

void EvtEtaPiPiGamma::init() 
{
	// check that there are 1 arguments
	checkNArg(1);
	checkNDaug(3);

	checkSpinParent(EvtSpinType::SCALAR);

	checkSpinDaughter(0,EvtSpinType::SCALAR);
	checkSpinDaughter(1,EvtSpinType::SCALAR);
	//checkSpinDaughter(2,EvtSpinType::VECTOR);
}


void EvtEtaPiPiGamma::initProbMax()
{

	setProbMax(2.2);

}

void EvtEtaPiPiGamma::decay(EvtParticle *p)  
{

	p->initializePhaseSpace(getNDaug(), getDaugs());

	EvtVector4R mompip = p->getDaug(0)->getP4();
	EvtVector4R mompim = p->getDaug(1)->getP4();
	//EvtVector4R momgamma = p->getDaug(2)->getP4();
	EvtVector4R mompipi = mompip + mompim;
	double masspip = p->getDaug(0)->mass();
	double masspim = p->getDaug(1)->mass();
	//double m_eta = p->mass();
	double spipi = mompipi.mass2();

	// extra kinematic phase space element factor
	double sigma3 = pow(1. - 4.*masspip*masspim/spipi, 3./2.);

	// additional decay factor
	double Fs = 1. + 2.12*spipi + 2.13*spipi*spipi + 13.80*spipi*spipi*spipi;

	// pull out parameters from arguments
	double alpha = getArg(0);

	EvtComplex amp(sqrt(sigma3 * Fs * (1. + alpha*spipi)),0.0);

	vertex(amp);

	return;

}


