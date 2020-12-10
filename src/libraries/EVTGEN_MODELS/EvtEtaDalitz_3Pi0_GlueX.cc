// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>

#include "EvtEtaDalitz_3Pi0_GlueX.h"

EvtEtaDalitz_3Pi0_GlueX::~EvtEtaDalitz_3Pi0_GlueX() {}

std::string EvtEtaDalitz_3Pi0_GlueX::getName(){

  return "ETA_DALITZ_3PI0_GLUEX";     

}


EvtDecayBase* EvtEtaDalitz_3Pi0_GlueX::clone(){

  return new EvtEtaDalitz_3Pi0_GlueX;

}

void EvtEtaDalitz_3Pi0_GlueX::init() 
{
	
	warning_counter =0;
	
	// check that there are 7 arguments
	checkNArg(3);
	checkNDaug(3);

	checkSpinParent(EvtSpinType::SCALAR);

	checkSpinDaughter(0,EvtSpinType::SCALAR);
	checkSpinDaughter(1,EvtSpinType::SCALAR);
	checkSpinDaughter(2,EvtSpinType::SCALAR);
}


void EvtEtaDalitz_3Pi0_GlueX::initProbMax()
{

	setProbMax(1.0);

}

void EvtEtaDalitz_3Pi0_GlueX::decay(EvtParticle *p)  
{

	p->initializePhaseSpace(getNDaug(), getDaugs());

	EvtVector4R mompi0_0 = p->getDaug(0)->getP4();
	EvtVector4R mompi0_1 = p->getDaug(1)->getP4();
	EvtVector4R mompi0_2 = p->getDaug(2)->getP4();
	double masspi0_0 = p->getDaug(0)->mass();
	double masspi0_1 = p->getDaug(1)->mass();
	double masspi0_2 = p->getDaug(2)->mass();
	double m_eta = p->mass();

	// Dalitz plot parameters (one recent reference: https://arxiv.org/abs/1803.02502)
	double Q = m_eta-masspi0_0-masspi0_1-masspi0_2;
	double x = sqrt(3.)*( (mompi0_1.get(0)-masspi0_1) - (mompi0_0.get(0)-masspi0_0) )/Q; 
	double y = (mompi0_2.get(0)-masspi0_2)*(3.0/Q)-1.0;
	double z = x*x + y*y;

	// pull out parameters from arguments
	double alpha = getArg(0); //Good measurements have been made for this parameter.
	double beta  = getArg(1); //Some ok measurements of this parameter.
	double gamma = getArg(2); //Barely probed with current experiments, but some values have been reported

	// Amplitude (square root of intensity/event distribution)
	EvtComplex amp(sqrt(1.0 + 2*alpha*z +2*beta*(3*x*x*y-y*y*y) +2*gamma*z*z ),0.0); // See Eq. 2 of reference

	// Jon experienced something weird where parameters sometimes fail to read in, but can't quite reproduce the problem.
	///// Add warning message, just in case this pops up again
	///// I assume there's really no reason to use all zero paramters, other than for debugging
	if( fabs(alpha)<0.00001 && fabs(beta)<0.00001  && fabs(gamma)<0.0001 && warning_counter < 1000) {
		std::cout << "Warning!!! Parameters are all zero. Results will be equivalent to phase space" << std::endl;
		std::cout << "If you think you supplied nonzero parameters, you might not be crazy. " << std::endl;
		std::cout << "Additional debugging is needed, but something about commented lines inside your eta particle decay portion of Decay.dec might be responsible" << std::endl;
		std::cout << std::endl;
		warning_counter++;
	}
		
	// std::cout << "Amp: " << abs2(amp) << std::endl;
	// std::cout << "z: " << z << std::endl;
	// std::cout << "alpha: " << alpha << std::endl;
	// std::cout << "beta: " << beta << std::endl;
	// std::cout << "gamma: " << gamma << std::endl << std::endl;

	vertex(amp);

	return;

}


