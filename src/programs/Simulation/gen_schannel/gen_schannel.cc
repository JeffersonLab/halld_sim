// Main program for generating eta events. 
#include "HDDM/hddm_s.h"
#include "particleType.h"

#include <TMath.h>
#include <TRandom3.h>
#include <TGenPhaseSpace.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <cstdlib>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm> 
#include <cctype>
#include <locale>
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

#include "UTILITIES/BeamProperties.h"

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// Masses
const double m_p=0.93827; // GeV
const double m_p_sq=m_p*m_p;

double zmin=50.0,zmax=80.0; // cm, target extent
int Nevents=10000;
int runNo=10000;
int generatedParticleType = 0;  // PDG type, default to an unknown type
bool debug=false;

// Diagnostic histograms
TH1D *cobrems_vs_E;

char output_file_name[250]="gen.hddm";
char beamConfigFile[250]="cobrems.conf";

void Usage(void){
  printf("gen_schannel: simple gamma + p -> X generator to be used as input for decay_evtgen.\n");
  printf(" Usage:  gen_schannel <options>\n");
  printf("   Options:  -N<number of events> (number of events to generate)\n");
  printf("             -B<beam.hddm>       (default: cobrems.conf)\n");
  printf("             -O<output.hddm>     (default: eta_gen.hddm)\n");
  printf("             -P<EvtGen particle number> (default: 0)\n");
  printf("             -R<run number>      (default: 10000)\n");
  printf("             -h                  (Print this message and exit.)\n");
  
  exit(0);
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
  int seed=0;
  if (narg==1){
    Usage();
  }
  for(int i=1; i<narg; i++){
    char *ptr = argv[i];
    
    if(ptr[0] == '-'){
      switch(ptr[1]){
      case 'h': Usage(); break;
      case 'O':
	sscanf(&ptr[2],"%s",output_file_name);
	break;
      case 'B':
	sscanf(&ptr[2],"%s",beamConfigFile);
	break;
      case 'P':
	sscanf(&ptr[2],"%d",&generatedParticleType);
	break;
      case 'N':
	sscanf(&ptr[2],"%d",&Nevents);
	break;
      case 'R':
	sscanf(&ptr[2],"%d",&runNo);
	break;
      case 'S':
	sscanf(&ptr[2],"%d",&seed);
	break;
      case 'd':
	debug=true;
	break;
      default:
	break;
      }
    }
  }
}

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
	// We don't actually care about this, just need the PDL...
 	EvtGen *myGenerator = new EvtGen(evtGenDecayFile.c_str(), evtGenParticleDefs.c_str(), eng,
  		     				 radCorrEngine, &extraModels);

	// open optional user decay file, if it exists
	struct stat buffer;   
	if(stat("userDecay.dec", &buffer) == 0)
	  	myGenerator->readUDecay("userDecay.dec");
 		     				 
}


// Put particle data into hddm format and output to file
void WriteEvent(unsigned int eventNumber, TLorentzVector &beam, float vert[3],
		int particle_type, TLorentzVector particle_vector, s_iostream_t *file)
{  
   s_PhysicsEvents_t* pes;
   s_Reactions_t* rs;
   s_Target_t* ta;
   s_Beam_t* be;
   s_Vertices_t* vs;
   s_Origin_t* origin;
   s_Products_t* ps;
   s_HDDM_t *thisOutputEvent = make_s_HDDM();
   thisOutputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
   pes->mult = 1;
   pes->in[0].runNo = runNo;
   pes->in[0].eventNo = eventNumber;
   pes->in[0].reactions = rs = make_s_Reactions(1);
   rs->mult = 1;
   // Beam 
   rs->in[0].beam = be = make_s_Beam();
   be->type = Gamma;
   be->properties = make_s_Properties();
   be->properties->charge = ParticleCharge(be->type);
   be->properties->mass = ParticleMass(be->type);
   be->momentum = make_s_Momentum();
   be->momentum->px = 0.;
   be->momentum->py = 0.;
   be->momentum->pz = beam.Pz();
   be->momentum->E  = beam.E();
   // Target
   rs->in[0].target = ta = make_s_Target();
   ta->type = Proton;
   ta->properties = make_s_Properties();
   ta->properties->charge = ParticleCharge(ta->type);
   ta->properties->mass = ParticleMass(ta->type);
   ta->momentum = make_s_Momentum();
   ta->momentum->px = 0.;
   ta->momentum->py = 0.;
   ta->momentum->pz = 0.;
   ta->momentum->E  = ParticleMass(ta->type);
   // Primary vertex 
   rs->in[0].vertices = vs = make_s_Vertices(1);
   vs->mult = 1;
   vs->in[0].origin = origin = make_s_Origin();
   vs->in[0].products = ps = make_s_Products(1);
   ps->mult = 1;
   origin->t = 0.0;
   origin->vx = vert[0];
   origin->vy = vert[1];
   origin->vz = vert[2];
   // Final state particle
   ps->in[0].type = UnknownParticle;  // zero out particle type info so that hdgeant won't decay the particle.  maybe there is a better way?
   //ps->in[0].type = particle_type;
   ps->in[0].pdgtype = particle_type;
   ps->in[0].id = 0; /* unique value for this particle within the event */
   ps->in[0].parentid = 0;  /* All internally generated particles have no parent */
   ps->in[0].mech = 0; // ???     
   ps->in[0].momentum = make_s_Momentum();
   ps->in[0].momentum->px = particle_vector.Px();
   ps->in[0].momentum->py = particle_vector.Py();
   ps->in[0].momentum->pz = particle_vector.Pz();
   ps->in[0].momentum->E  = particle_vector.E();

   flush_s_HDDM(thisOutputEvent,file);
}

// Create some diagnostic histograms
void CreateHistograms(string beamConfigFile)
{
  BeamProperties beamProp(beamConfigFile);
  cobrems_vs_E = (TH1D*)beamProp.GetFlux();
}


//-----------
// main
//-----------
int main(int narg, char *argv[])
{  
  ParseCommandLineArguments(narg, argv);

  InitEvtGen();

  // open HDDM file
  s_iostream_t *file = init_s_HDDM(output_file_name);
 
  // Initialize random number generator
  //TRandom3 *myrand=new TRandom3(0);// If seed is 0, the seed is automatically computed via a TUUID object, according to TRandom3 documentation

  // Fixed target
  TLorentzVector target(0.,0.,0.,m_p);

  //----------------------------------------------------------------------------
  // Get production (Egamma range) and decay parameters from input file
  //----------------------------------------------------------------------------

  // Get beam properties configuration file
  cout << "Photon beam configuration file " << beamConfigFile << endl;

  // Create some diagonistic histographs
  CreateHistograms(beamConfigFile);

  cout << "Generating events..." << endl;

  //----------------------------------------------------------------------------
  // Event generation loop
  //----------------------------------------------------------------------------
  for (int i=1;i<=Nevents;i++) {
    double Egamma=0.;

    // vertex position at target
    float vert[4]={0.,0.,0.,0.};

    // use the rejection method to produce eta's based on the cross section
	// First generate a beam photon using bremsstrahlung spectrum
	Egamma = cobrems_vs_E->GetRandom();

	// CM energy
	double s=m_p*(m_p+2.*Egamma);
	double Ecm=sqrt(s);

    // beam 4-vector (ignoring px and py, which are extremely small)
    TLorentzVector beam(0.,0.,Egamma,Egamma);

	// generated particle
	double mass = EvtPDL::getMass(EvtPDL::evtIdFromStdHep(generatedParticleType));
	double Pz = sqrt(Ecm*Ecm - mass*mass);
	TLorentzVector part4v(0.,0.,Pz,Ecm);

    // Write Event to HDDM file
    WriteEvent(i,beam,vert,generatedParticleType,part4v,file);
    
    //if ((i%(Nevents/10))==0) cout << 100.*double(i)/double(Nevents) << "\% done" << endl;
  }


  // Close HDDM file
  close_s_HDDM(file);
  cout<<endl<<"Closed HDDM file"<<endl;
  cout<<" "<<Nevents<<" event written to "<<output_file_name<<endl;


  return 0;
}
