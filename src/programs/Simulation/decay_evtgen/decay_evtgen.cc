// decay_evtgen.cc
// Description
// Sean Dobbs, sdobbs@fsu.edu (2019)

#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

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

#include "particleType.h"
#include "HDDM/hddm_s.hpp"
#include "EVTGEN_MODELS/RegisterGlueXModels.h"

#include "TLorentzVector.h"
#include "TVector3.h"

typedef struct {
	int id = -1;
	int type = -1;   // PDG type
	bool decayed = false;
	TLorentzVector momentum;
	TVector3 vertex;
} gen_particle_info_t;

string INPUT_FILE = "";
string OUTPUT_FILE = "";
EvtGen *myGenerator = nullptr;

void InitEvtGen();
void ParseCommandLineArguments(int narg,char *argv[]);
void Usage(void);
void ParseVertices(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info);
void DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info);

//-------------------------------
// InitEvtGen
//-------------------------------
void InitEvtGen()
{
  	// initialize EvtGen
  	const char* evtgen_home_env_ptr = std::getenv("EVTGEN_HOME");
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
	if(stat("userDecay.dec", &buffer) == 0)
	  	myGenerator->readUDecay("userDecay.dec");


}


//-------------------------------
// ParseVertices
//-------------------------------
void ParseVertices(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, int &max_particle_id)
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
			part_info.type = it_product->getPdgtype();
			part_info.decayed = false;
			TLorentzVector mom(it_product->getMomentum().getPx(), it_product->getMomentum().getPy(), 
								it_product->getMomentum().getPz(), it_product->getMomentum().getE());
			part_info.momentum = mom;
			TVector3 vertex(it_vertex->getOrigin().getVx(), it_vertex->getOrigin().getVy(),
							it_vertex->getOrigin().getVz());
			part_info.vertex = vertex;
			
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
void DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, 
					int &max_particle_id)
{
	EvtParticle* parent(0);
	for(auto &part : particle_info) {
		if(part.decayed)
			continue;
	
		// Set up the parent particle
		EvtId partId = EvtPDL::getId(std::string(EvtGenString(PDGtoPType(part.type))));
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
		switch(PDGtoPType(part.type)) {
			case KShort: case KLong: case Lambda: case SigmaPlus:
  			case Sigma0: case SigmaMinus: case Xi0: case XiMinus:  case OmegaMinus:
  			case AntiLambda: case AntiSigmaMinus: case AntiSigma0: case AntiSigmaPlus:
  			case AntiXi0: case AntiXiPlus: case AntiOmegaPlus:
  				continue;
			default:
				myGenerator->generateDecay(parent);
		}
		
		// add decay vertex and daughter particles to HDDM record, if a decay happened
		if(parent->getNDaug() > 0) {
			part.decayed = true;    
			hddm_s::ReactionList rs = hddmevent->getReactions();
			hddm_s::VertexList vertices = rs().addVertices();
			hddm_s::ProductList ps = vertices().addProducts(parent->getNDaug());
			
			// set the event vertex
		    hddm_s::OriginList os = vertices().addOrigins();
			os().setT(0.0);
			os().setVx(part.vertex.x());
			os().setVy(part.vertex.y());
			os().setVz(part.vertex.z());
			
			// save the information on the daughter particles
			vector< gen_particle_info_t > decay_particle_info;
			for (unsigned int i=0; i<parent->getNDaug(); i++){
				ps(i).setType(PDGtoPType(parent->getDaug(i)->getPDGId()));
				ps(i).setPdgtype(parent->getDaug(i)->getPDGId());
				ps(i).setId(++max_particle_id);   // unique value for this particle within the event
				ps(i).setParentid(part.id);       // set ID of parent particle
				ps(i).setMech(0);        // ???     
		        hddm_s::MomentumList pmoms = ps(i).addMomenta();
				pmoms().setPx(parent->getDaug(i)->getP4Lab().get(1));
				pmoms().setPy(parent->getDaug(i)->getP4Lab().get(2));
				pmoms().setPz(parent->getDaug(i)->getP4Lab().get(3));
				pmoms().setE(parent->getDaug(i)->getP4Lab().get(0));
				
				// save the same info so that we can go through these particles and see if they need to decay
				gen_particle_info_t part_info;
				part_info.id = max_particle_id;
				part_info.type = parent->getDaug(i)->getPDGId();
				part_info.decayed = false;
				TLorentzVector mom(parent->getDaug(i)->getP4Lab().get(1), parent->getDaug(i)->getP4Lab().get(2), 
								   parent->getDaug(i)->getP4Lab().get(3), parent->getDaug(i)->getP4Lab().get(0));
				part_info.momentum = mom;
				part_info.vertex = part.vertex;
			
				decay_particle_info.push_back(part_info);    
			}
		
			// go ahead and decay the particles we just generated
			DecayParticles(hddmevent, decay_particle_info, max_particle_id);
		}
		
		parent->deleteTree();   //cleanup
	}

}


//-------------------------------
// main
//-------------------------------ƒ
int main(int narg, char *argv[])
{
	ParseCommandLineArguments(narg,argv);

	if (INPUT_FILE == "") {
	  cerr << "No input file!" << endl;
	}

	// Open input file
	ifstream *infile = new ifstream(INPUT_FILE);
	if (! infile->is_open()) {
	  cerr << "Unable to open file \"" << INPUT_FILE << "\" for reading."
				<< endl;
	  exit(-2);
	}
	cout << "Opening Input File:  " << INPUT_FILE << " ..." << endl;
	hddm_s::istream *instream = new hddm_s::istream(*infile);

	// Open output file
	ofstream *outfile = new ofstream(OUTPUT_FILE.c_str());
	if (! outfile->is_open()) {
	  cerr << "Unable to open output file \"" << OUTPUT_FILE
				<< "\" for writing." << endl;
	  exit(-3);
	}
	cout << "Opening Output File:  " << OUTPUT_FILE << " ..." << endl;
	hddm_s::ostream *outstream = new hddm_s::ostream(*outfile);

	InitEvtGen();
	
	int event_count = 1;
	hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
	while(*instream >> *hddmevent) {
		int num_particles = -1;   // number of particles in the event
		int max_particle_id = 0;  // needed for generating decay particles

		if( (event_count++%1000) == 0) {
			cout << "Processed " << event_count << " events ..." << endl;
		}

		vector< gen_particle_info_t > particle_info;
		ParseVertices(hddmevent, particle_info, max_particle_id);    // fill particle info vector
		DecayParticles(hddmevent, particle_info, max_particle_id);   // run EvtGen decays based on particle info vector
	
	   	*outstream << *hddmevent;  // save event
	}

	// cleanup
	delete instream;
	delete outstream;

	return 0;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
   if (narg < 2) {
      Usage();
      exit(0);
   }

   for(int i=1; i < narg; i++) {
      if (argv[i][0]=='-') {
         char *ptr = &argv[i][1];
         switch(*ptr) {
            case 'o':
              OUTPUT_FILE = &ptr[1];
              break;
            default:
              cerr << "Unknown option \"" << argv[i] << "\"" << endl;
              Usage();
              exit(-1);
         }
      }
      else {
         INPUT_FILE = argv[i];
      }
   }
   
   if(OUTPUT_FILE == "") {
	   // Determine output filename from input filename
	   OUTPUT_FILE = INPUT_FILE;
	   size_t pos = OUTPUT_FILE.find_last_of(".");
	   if (pos != string::npos) OUTPUT_FILE.erase(pos);
	   OUTPUT_FILE += "_decayed.hddm";
   }
}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
  cout << endl;
  cout << "Usage:" << endl;
  cout << "       decay_evtgen [options] file.hddm" << endl;
  cout << endl;
  cout << "Decay any particles via EvtGen (update!)" << endl;
  cout << endl;
  cout << " options:" << endl;
  cout << endl;
  cout << "  -o\"output_file_name\"    "
               "set the file name used for output." << endl;
  cout << "  -h                        "
               "print this usage statement." << endl;
  cout << endl;
}

