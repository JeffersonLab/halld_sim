
#include <EvtGenDecayer.h>

#include <particleType.h>

#include <sys/stat.h>

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

#include "EVTGEN_MODELS/RegisterGlueXModels.h"




////////////////////////////////////////////////////
// EvtGenDecayer
//
// Initialization (constructor)
////////////////////////////////////////////////////
EvtGenDecayer::EvtGenDecayer()
{
	// TODO: Make it so that some of these options can be passed in as arguments
    string USER_DECAY = "userDecay.dec";

  	// initialize EvtGen
  	const char* evtgen_home_env_ptr = std::getenv("EVTGENDIR");
  	string EVTGEN_HOME = (evtgen_home_env_ptr==nullptr) ? "." : evtgen_home_env_ptr;  // default to the current directory
  	
    // Define the random number generator
    EvtRandomEngine* eng = 0;
#ifdef EVTGEN_CPP11
 	 // Use the Mersenne-Twister generator (C++11 only)
 	 eng = new EvtMTRandomEngine();
#else
 	 eng = new EvtSimpleRandomEngine();
#endif

 	 EvtRandom::setRandomEngine(eng);

 	 EvtAbsRadCorr* radCorrEngine = 0;
 	 std::list<EvtDecayBase*> extraModels;

#ifdef EVTGEN_EXTERNAL
     // this is all required for using PHOTOS
 	 bool convertPythiaCodes(false);
	 bool useEvtGenRandom(true);
 	 EvtExternalGenList genList(convertPythiaCodes, "", "gamma", useEvtGenRandom);
	 radCorrEngine = genList.getPhotosModel();
	 extraModels = genList.getListOfModels();
#endif

 	//Initialize the generator - read in the decay table and particle properties
  	const char* evtgen_decay_file_ptr = std::getenv("EVTGEN_DECAY_FILE");
 	string evtGenDecayFile = (evtgen_decay_file_ptr==nullptr) ? EVTGEN_HOME + "/DECAY.DEC" : evtgen_decay_file_ptr;
  	const char* evtgen_particle_defs_ptr = std::getenv("EVTGEN_PARTICLE_DEFINITIONS");
 	string evtGenParticleDefs = (evtgen_particle_defs_ptr==nullptr) ? EVTGEN_HOME + "/evt.pdl" : evtgen_particle_defs_ptr;
 	myGenerator = new EvtGen(evtGenDecayFile.c_str(), evtGenParticleDefs.c_str(), eng,
  		     				 radCorrEngine, &extraModels);
  		     				 
  	// Make special GlueX definitions
  	GlueX_EvtGen::RegisterGlueXModels();

	// open optional user decay file, if it exists
	struct stat buffer;   
	if(stat(USER_DECAY.c_str(), &buffer) == 0)
	  	myGenerator->readUDecay(USER_DECAY.c_str());


}

////////////////////////////////////////////////////
// decayParticle
//
// inputs: 4-vector of parent particle and particleID (in Particle_t format)
// output: vector of (TLorentzVector, Particle_t) pairs
////////////////////////////////////////////////////
vector< pair< TLorentzVector, int > > EvtGenDecayer::decayParticle( const TLorentzVector& parent4v, int particleID )
{
	// Set up the parent particle in EvtGen format
	EvtParticle* parent(0);

	// this will work as long as we are dealing with particles with a well-defined Particle_t value
	// could be more flexible to handle more PDG types...
	EvtId partId = EvtPDL::evtIdFromStdHep(PDGtype(static_cast<Particle_t>(particleID)));

//  is this functionality that we want?					
// 	// allow an optional remapping of particle types based on the PDG type, to deal with particles
// 	// that are not properly defined in the GlueX framework
// 	int current_pdgtype = EvtPDL::getStdHep( partId );
// 	if(PDGTYPE_CONVERSION_MAP.find(current_pdgtype) != PDGTYPE_CONVERSION_MAP.end()) {
// 		//int oldid = current_pdgtype;
// 		partId = EvtPDL::evtIdFromStdHep(PDGTYPE_CONVERSION_MAP[current_pdgtype]);
// 		//int newid = EvtPDL::getStdHep( partId );
// 	}
	
	EvtVector4R pInit(parent4v.E(), parent4v.Px(), parent4v.Py(), parent4v.Pz());
	parent = EvtParticleFactory::particleFactory(partId, pInit);

	// Generate the particle decays - this is where the real meat happens
	// This generates the particle decay based on the predefined/user-defined decays
	// If the particle is "stable", then there is no defined decay in the default files
	// which is true for particles like the electron and proton, but also the pion and kaon.
	// "long-lived" particles should really be decayed by Geant.  
	// I need to add an exclusion list so that we don't decay neutral kaons, hyperons, etc.
	// hardcode it for now...
	switch(particleID) {
		case KShort: case KLong: case Lambda: case SigmaPlus:
		case Sigma0: case SigmaMinus: case Xi0: case XiMinus:  case OmegaMinus:
		case AntiLambda: case AntiSigmaMinus: case AntiSigma0: case AntiSigmaPlus:
		case AntiXi0: case AntiXiPlus: case AntiOmegaPlus:
			;  // no-op
		default:
			myGenerator->generateDecay(parent);
	}

	// build the output
	vector< pair< TLorentzVector, int > > children;
	
	if(parent->getNDaug() > 0) {
		for (unsigned int i=0; i<parent->getNDaug(); i++) {
			TLorentzVector part4v(parent->getDaug(i)->getP4Lab().get(1), parent->getDaug(i)->getP4Lab().get(2), 
					parent->getDaug(i)->getP4Lab().get(3), parent->getDaug(i)->getP4Lab().get(0));
			children.push_back(pair< TLorentzVector, int >(part4v, PDGtoPType(parent->getDaug(i)->getPDGId()) ));
		}
	}
	
	parent->deleteTree();   //cleanup
	
	return children;

}
