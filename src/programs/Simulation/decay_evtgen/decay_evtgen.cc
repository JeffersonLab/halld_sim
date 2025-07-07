// decay_evtgen.cc
// Description
// Sean Dobbs, sdobbs@fsu.edu (2019)

#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
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

//#include "particleType.h"  // manage this function locally to cover older halld_recon versions
#include "HDDM/hddm_s.hpp"
#include "EVTGEN_MODELS/RegisterGlueXModels.h"
#include "evtgenParticleString.h"

#include "TLorentzVector.h"
#include "TVector3.h"

typedef struct {
	int id = -1;
	Particle_t type = UnknownParticle;
	int pdgtype = -1;   
	bool decayed = false;
	TLorentzVector momentum;
	TVector3 vertex;
	hddm_s::ProductList::iterator hddmProduct;
} gen_particle_info_t;

string INPUT_FILE = "";
string OUTPUT_FILE = "";
string USER_DECAY = "userDecay.dec";
EvtGen *myGenerator = nullptr;

bool PROCESS_ALL_EVENTS = true;
int NUM_EVENTS_TO_PROCESS = -1;
bool GEN_SCHANNEL = false;

map<int,int> PDGTYPE_CONVERSION_MAP;

void InitEvtGen();
void ParseCommandLineArguments(int narg,char *argv[]);
void Usage(void);
void ParseVertices(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info);
void DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, int &vertex_id);

//-------------------------------
// InitEvtGen
//-------------------------------
void InitEvtGen()
{
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
void DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, 
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
					part.hddmProduct->setType(UnknownParticle);  // zero out particle type info so that hdgeant won't decay the particle.  maybe there is a better way?
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
			os().setT(0.0);
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
	  // next line commented out, an unused variable
	  //		int num_particles = -1;   // number of particles in the event
		int max_particle_id = 0;  // needed for generating decay particles

		if( (event_count++%1000) == 0) {
			cout << "Processed " << event_count << " events ..." << endl;
		}

		vector< gen_particle_info_t > particle_info;
		int vertex_id = 0;
		ParseVertices(hddmevent, particle_info, max_particle_id);    // fill particle info vector
		DecayParticles(hddmevent, particle_info, max_particle_id, vertex_id);   // run EvtGen decays based on particle info vector
	
	   	*outstream << *hddmevent;  // save event

		// see if we should stop processing
		if(!PROCESS_ALL_EVENTS) {
		  if(event_count >= NUM_EVENTS_TO_PROCESS)
		    break;
		}
	}

	// cleanup
	delete instream;
	delete outstream;

	return 0;
}

//-------------------------------
// ConvertStringInt
//-------------------------------
bool ConvertStringInt(string &the_str, int &out_int) 
{
	istringstream the_istream(the_str);
	the_istream >> out_int;
	if(the_istream.fail())
		return false;

	return true;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[])
{
  string num_events_str;
  size_t seperator_index;
  string argstr;
  int pdgtype1, pdgtype2;
  string substr1, substr2;

   if (narg < 2) {
      Usage();
      exit(0);
   }

   for(int i=1; i < narg; i++) {
      if (argv[i][0]=='-') {
         char *ptr = &argv[i][1];
         switch(*ptr) {
            case 'n':
	      PROCESS_ALL_EVENTS = false;
	      num_events_str = &ptr[1];
	      NUM_EVENTS_TO_PROCESS = std::stoi(num_events_str);
              break;
            case 'o':
              OUTPUT_FILE = &ptr[1];
              break;
            case 'u':
              USER_DECAY = &ptr[1];
              break;
            case 'S':
              GEN_SCHANNEL = true;
              break;
            case 'X':
              // assume input of the form -XN_M
              // where N and M are both PDG ID numbers
              argstr = &ptr[1];
              seperator_index = argstr.find("_");
              if(seperator_index == string::npos) {
              	cerr << " Invalid -X format: " << argstr << endl;
              } else {
              	// let's actually try parsing
              	substr1 = argstr.substr(0,seperator_index);
              	if(!ConvertStringInt(substr1,pdgtype1)) {
              		cerr << " Invalid -X format: " << argstr << " - bad type " << substr1 << endl;
              	} else {
              		substr2 = argstr.substr(seperator_index+1);
              		if(!ConvertStringInt(substr2,pdgtype2)) {
              	 		cerr << " Invalid -X format: " << argstr << " - bad type " << substr2 << endl;
              	 	} else {
              	 		cerr << "setting up conversion of type " << pdgtype1 << " to " << pdgtype2 << endl;
						PDGTYPE_CONVERSION_MAP[pdgtype1] = pdgtype2;          	 
              	 	}
              	 }
              }
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
  cout << "Decay thrown particles via EvtGen" << endl;
  cout << "The particle types and decays are defined in $EVTGENDIR/evt.pdl " << endl
       << "  and $EVTGENDIR/DECAY.DEC respectively." << endl;
  cout << "The location of these files can be override by the environment variables " << endl
       << "  EVTGEN_PARTICLE_DEFINITIONS and EVTGEN_DECAY_FILE respectively." << endl;
  cout << endl << "For more documentation on EvtGen, go to https://evtgen.hepforge.org/" << endl;
  cout << endl;
  cout << " options:" << endl;
  cout << endl;
  cout << "  -nNumEvents               "
               "number of events to process (default: all)" << endl;
  cout << "  -o\"output_file_name\"    "
               "set the file name used for output (default: append \"_decayed\")" << endl;
  cout << "  -u\"user_decay_file_name\"    "
               "set the file name of the user decay file (default: userDecay.dec)" << endl;
  cout << "  -XM_N                         "
      "translate particles of PDG ID M to those of ID N.  This option can be specified multiple times" << endl;
  cout << "  -h                        "
               "print this usage statement." << endl;
  cout << endl;
}

