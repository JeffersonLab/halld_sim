#include "EvtGenInterface.h"


//-------------------------------
// InitEvtGen
//-------------------------------
void EvtGenInterface::InitEvtGen()
{
	// Do some sanity checks here so that we don't initialize twice!!

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


//-------------------------------
// Decay
//-------------------------------
void EvtGenInterface::Decay(hddm_s::HDDM * hddmevent)
{
	int max_particle_id = 0;  // needed for generating decay particles
	vector< gen_particle_info_t > particle_info;
	int vertex_id = 0;
	
	ParseVertices(hddmevent, particle_info, max_particle_id);    // fill particle info vector
	DecayParticles(hddmevent, particle_info, max_particle_id, vertex_id);   // run EvtGen decays based on particle info vector
}



//-------------------------------
// ParseVertices
//-------------------------------
void EvtGenInterface::ParseVertices(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, int &max_particle_id)
{
   hddm_s::VertexList vertices = hddmevent->getVertices();
   hddm_s::VertexList::iterator it_vertex;

	// put the thrown particles in the HDDM record into some intermediate format
	for (it_vertex = vertices.begin(); it_vertex != vertices.end(); ++it_vertex) {
    	hddm_s::ProductList &products = it_vertex->getProducts();
      	hddm_s::ProductList::iterator it_product;
      	for (it_product = products.begin(); it_product != products.end(); ++it_product) {
        	// ignore intermediaries in the MC record
			//if (it_product->getType() <= 0)
         	//continue;

			gen_particle_info_t part_info;
			part_info.id = it_product->getId();
			part_info.type = it_product->getType();
			part_info.pdgtype = it_product->getPdgtype();
			part_info.decayed = false;
			TLorentzVector mom(it_product->getMomentum().getPx(), it_product->getMomentum().getPy(), 
								it_product->getMomentum().getPz(), it_product->getMomentum().getE());
			part_info.momentum = mom;
			TVector3 vertex(it_vertex->getOrigin().getVx(), it_vertex->getOrigin().getVy(),
							it_vertex->getOrigin().getVz());
			part_info.vertex = vertex;
			part_info.hddmProduct = it_product;
			
			// track the maximum particle ID for when we need to add more particles
			if(max_particle_id < part_info.id)
				max_particle_id = part_info.id;
			
			// check parent and flag that is has been decayed
			int parent_id = it_product->getParentid();
			if(parent_id) {
				for(auto &part : particle_info) {
					if(part.id == parent_id) {
						part.decayed = true;
						break;
					}
				}
			}
			
			particle_info.push_back(part_info);    
		}
    }

}

//-------------------------------
// DecayParticles
//-------------------------------
void EvtGenInterface::DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, 
					int &max_particle_id, int &vertex_id)
{
	EvtParticle* parent(0);
	for(auto &part : particle_info) {
//   		cerr << "part id = " << part.id << "  type = " << part.type 
//  			<< "  pdgtype = " << part.pdgtype << "  evtid = " <<   EvtPDL::evtIdFromStdHep(part.pdgtype) << endl;
		
		if(part.decayed)
			continue;
		if(!GEN_SCHANNEL && part.type == 0) {
			cout << "Particle of type 0 detected, skipping!" << endl;
			continue;
		}
	
		// Set up the parent particle
		//EvtId partId = EvtPDL::getId(std::string(EvtGenString(part.type)));  // do we really want to use EvtGenString?  maybe use some internal EvtGen function instead...
		EvtId partId;   
		//  do some mapping between the particle IDs in the HDDM files and EvtGen IDs
		//  this gets complicated since most generators are good about writing out particle types
		//  in the PDG scheme, but some aren't.  should fix that at some point
		if(part.pdgtype != 0)
			// prefer determining particle type from its PDG numbering
			partId = EvtPDL::evtIdFromStdHep(part.pdgtype);  
		else
			// use the standard particle type as a backup - would be nice
			// if we didn't have these dependencies!
			// note that bggen-derived generators tend to have this issue
			partId = EvtPDL::getId(std::string(EvtGenOutputString(part.type))); 
						
		// allow an optional remapping of particle types based on the PDG type, to deal with particles
		// that are not properly defined in the GlueX framework
		int current_pdgtype = EvtPDL::getStdHep( partId );
		if(PDGTYPE_CONVERSION_MAP.find(current_pdgtype) != PDGTYPE_CONVERSION_MAP.end()) {
			//int oldid = current_pdgtype;
			partId = EvtPDL::evtIdFromStdHep(PDGTYPE_CONVERSION_MAP[current_pdgtype]);
			//int newid = EvtPDL::getStdHep( partId );
		}
		
		EvtVector4R pInit(part.momentum.E(), part.momentum.Px(), 
							part.momentum.Py(), part.momentum.Pz());
		parent = EvtParticleFactory::particleFactory(partId, pInit);

		// Generate the particle decays - this is where the real meat happens
		// This generates the particle decay based on the predefined/user-defined decays
		// If the particle is "stable", then there is no defined decay in the default files
		// which is true for particles like the electron and proton, but also the pion and kaon.
		// "long-lived" particles should really be decayed by Geant.  
		// I need to add an exclusion list so that we don't decay neutral kaons, hyperons, etc.
		// hardcode it for now...
		switch(part.type) {
			case KShort: case KLong: case Lambda: case SigmaPlus:
  			case Sigma0: case SigmaMinus: case Xi0: case XiMinus:  case OmegaMinus:
  			case AntiLambda: case AntiSigmaMinus: case AntiSigma0: case AntiSigmaPlus:
  			case AntiXi0: case AntiXiPlus: case AntiOmegaPlus:
  				continue;
			default:
				myGenerator->generateDecay(parent);
				if(parent->getNDaug() > 0) 
					part.hddmProduct->setType(Unknown);  // zero out particle type info so that hdgeant won't decay the particle.  maybe there is a better way?
				break;
		}
		
		// add decay vertex and daughter particles to HDDM record, if a decay happened
		if(parent->getNDaug() > 0) {
			part.decayed = true;    
			hddm_s::ReactionList rs = hddmevent->getReactions();
			hddm_s::VertexList vertices = rs().addVertices();
			hddm_s::ProductList ps = vertices().addProducts(parent->getNDaug());
			hddm_s::ProductList::iterator it_product = ps.begin();

			// set the event vertex
		    hddm_s::OriginList os = vertices().addOrigins();
			os().setT(0.0);  // probably need to correctly set this?
			os().setVx(part.vertex.x());
			os().setVy(part.vertex.y());
			os().setVz(part.vertex.z());
			vertex_id++;
			
			// save the information on the daughter particles
			vector< gen_particle_info_t > decay_particle_info;
			for (unsigned int i=0; i<parent->getNDaug(); i++, it_product++){
				ps(i).setDecayVertex(vertex_id);
				ps(i).setType(PDGtoPType(parent->getDaug(i)->getPDGId()));
				ps(i).setPdgtype(parent->getDaug(i)->getPDGId());
				ps(i).setId(++max_particle_id);   // unique value for this particle within the event
				ps(i).setParentid(part.id);       // set ID of parent particle
				ps(i).setMech(0);        // ???     
				hddm_s::MomentumList pmoms = ps(i).addMomenta();
				//pmoms().setPx(parent->getDaug(i)->getP4Lab().get(1));
				//pmoms().setPy(parent->getDaug(i)->getP4Lab().get(2));
				//pmoms().setPz(parent->getDaug(i)->getP4Lab().get(3));
				// recalculate the energy to make sure we get the "correct" mass without round-off errors, since
				// HDDM stores floating point numbers as floats
				float px = parent->getDaug(i)->getP4Lab().get(1);
				float py = parent->getDaug(i)->getP4Lab().get(2);
				float pz = parent->getDaug(i)->getP4Lab().get(3);
				float mass =  EvtPDL::getMass(EvtPDL::evtIdFromStdHep(parent->getDaug(i)->getPDGId()));
				float E = sqrt(mass*mass + px*px + py*py + pz*pz);
				pmoms().setPx(px);
				pmoms().setPy(py);
				pmoms().setPz(pz);
				pmoms().setE(E);
				
				// save the same info so that we can go through these particles and see if they need to decay
				gen_particle_info_t part_info;
				part_info.id = max_particle_id;
				part_info.type = PDGtoPType(parent->getDaug(i)->getPDGId());
				part_info.pdgtype = parent->getDaug(i)->getPDGId();
				part_info.decayed = false;
				TLorentzVector mom(parent->getDaug(i)->getP4Lab().get(1), parent->getDaug(i)->getP4Lab().get(2), 
								   parent->getDaug(i)->getP4Lab().get(3), parent->getDaug(i)->getP4Lab().get(0));
				part_info.momentum = mom;
				part_info.vertex = part.vertex;
				part_info.hddmProduct = it_product;
			
				decay_particle_info.push_back(part_info);    
			}
		
			// go ahead and decay the particles we just generated
			DecayParticles(hddmevent, decay_particle_info, max_particle_id, vertex_id);
		}
		
		parent->deleteTree();   //cleanup
	}

}


