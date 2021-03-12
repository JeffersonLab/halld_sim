#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <math.h>
#include "HDDM/hddm_s.hpp"
#include "particleType.h"
using namespace std;

void ParseCommandLineArguments(int narg,char *argv[]);
void open_ascii_input(void);
void open_hddm_output(void);
void close_files(void);
void clean_arr(void);
int Str2GeantParticleID(char *str);
void Usage(void);
void parse_kinematic_values_from_ASCII(void);
void get_topology(void);
void get_vt(void);
void copy_arr_to_product(hddm_s::ProductList ps, int ps_id, int loc_idx);
void write_hddm_events(int runNumber, int eventNumber);
bool final_states_flag(int loc_idx);


char *INPUT_FILE=NULL;
string OUTPUT_FILE("output.hddm");
std::ofstream *outfile = NULL;
hddm_s::ostream *outstream = NULL;
ifstream *infile=NULL;

Particle_t targetType = Proton;
Particle_t beamType = Gamma;
int USER_RUNNUMBER = 10000;

//array to store values from ASCII 
int nParticles, mechFlag;
float beamEnergy;
std::vector<int> arr_Particletype;
std::vector<int> arr_id, arr_parentid;
std::vector<float> arr_E, arr_px, arr_py, arr_pz, arr_x, arr_y, arr_z, arr_vx, arr_vy, arr_vz, arr_vt, arr_t;

// containers for processed information
std::map<int, set<int>> loc_topology_map;

//-------------------------------
// main
//-------------------------------
int main(int narg, char *argv[])
{
  ParseCommandLineArguments(narg,argv);
   
  // Open input file
  open_ascii_input();

  //open hddm file to write
  open_hddm_output();

  // Loop over events
  int eventNumber = 0;
   
  while (! infile->eof()) 
  {
      //clean up memory from last event
      clean_arr();
      nParticles=0;
      beamEnergy = 0.0;
      mechFlag = 0;
      eventNumber++;

      //set up run numbers
      int runNumber = 10000; //DUMMY
      if(USER_RUNNUMBER != 0) runNumber = USER_RUNNUMBER;

      // read kinematic values from the ascii
      parse_kinematic_values_from_ASCII();

      // get the topology
      get_topology();

      //get flytime
      get_vt();

      // write in the hddm format
      write_hddm_events(runNumber, eventNumber);

    }//end of the while loop over all events

   // Close output file
   close_files();
   std::cout << "Wrote " << eventNumber << " events to " << OUTPUT_FILE
            << std::endl;
   
   return 0;
}


   








//------------------
// write_hddm_events
//------------------
void write_hddm_events(int runNumber, int eventNumber){
  // Start a new event
  hddm_s::HDDM record;
  hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
  pes().setRunNo(runNumber);
  pes().setEventNo(eventNumber);

  // Reaction
  hddm_s::ReactionList rs = pes().addReactions();

  // Target
  hddm_s::TargetList ts = rs().addTargets();
  ts().setType(targetType);
  hddm_s::PropertiesList tpros = ts().addPropertiesList();
  tpros().setCharge(ParticleCharge(targetType));
  tpros().setMass(ParticleMass(targetType));
  hddm_s::MomentumList tmoms = ts().addMomenta();
  tmoms().setPx(0);
  tmoms().setPy(0);
  tmoms().setPz(0);
  tmoms().setE(ParticleMass(targetType));

  // Beam
  hddm_s::BeamList bs = rs().addBeams();
  bs().setType(beamType);
  hddm_s::PropertiesList bpros = bs().addPropertiesList();
  bpros().setCharge(ParticleCharge(beamType));
  bpros().setMass(ParticleMass(beamType));
  hddm_s::MomentumList bmoms = bs().addMomenta();
  bmoms().setPx(-tmoms().getPx());
  bmoms().setPy(-tmoms().getPy());
  bmoms().setPz(beamEnergy);
  bmoms().setE(beamEnergy); 

  // create list of vertices
  int n_vertices = (int)loc_topology_map.size();
  hddm_s::VertexList vs = rs().addVertices(n_vertices); 

  //loop over the topology map to set up the vertex
  int loc_vertex_idx = 0;
  for(auto it = loc_topology_map.begin(); it != loc_topology_map.end(); it++){

      hddm_s::OriginList loc_os = vs(loc_vertex_idx).addOrigins();

      // create list of products
      int nParts = (int)(it->second.size());
      hddm_s::ProductList ps = vs(loc_vertex_idx).addProducts(nParts);

      //set up particles in the ProductList
      int ps_id = 0;
      for(auto itt = it->second.begin(); itt != it->second.end(); ++itt){
        int loc_particle_id = *itt;
        copy_arr_to_product(ps, ps_id, loc_particle_id);
        ps_id++;
      }

      //set up the vertex position
      int loc_particle_id = (int)(*it->second.begin());
      loc_os().setT(arr_vt.at(loc_particle_id-1));
      loc_os().setVx(arr_x.at(loc_particle_id-1));
      loc_os().setVy(arr_y.at(loc_particle_id-1));
      loc_os().setVz(arr_z.at(loc_particle_id-1));
      
      


      loc_vertex_idx++;
  }

  //display number of events processed along the way
  if (nParticles > 0){
     *outstream << record;
     if (eventNumber%10000 == 0) {
        std::cout << "Wrote event " << eventNumber << "\r";
        std::cout.flush();
     }
  }

}

                           
//-----------------
// copy_arr_to_product
//-----------------
void copy_arr_to_product(hddm_s::ProductList ps, int ps_id, int loc_idx){
    // check for intermediate particle (i.e. ones that GEANT does not track)
    // set their "type" explicitly to 0 to tell hdgeant to ignore them. 
    int type;
    if(!final_states_flag(loc_idx)){
       type = 0;
    }
    else{
      type = arr_Particletype.at(loc_idx-1);
    }
    // write in the product profile
    ps(ps_id).setType((Particle_t)type);
    ps(ps_id).setPdgtype(PDGtype((Particle_t)type));
    ps(ps_id).setId(loc_idx);         // unique value for this particle within the event //
    ps(ps_id).setParentid(arr_parentid.at(loc_idx-1));     // All internally generated particles have no parent //
    ps(ps_id).setMech(mechFlag);         // maybe this should be set to something? //
    hddm_s::MomentumList pmoms = ps(ps_id).addMomenta();  
                         pmoms().setPx(arr_px.at(loc_idx-1));
                         pmoms().setPy(arr_py.at(loc_idx-1));
                         pmoms().setPz(arr_pz.at(loc_idx-1));
                         pmoms().setE(arr_E.at(loc_idx-1));            

}

bool final_states_flag(int loc_idx){

  if((int)loc_topology_map.size() == 1){
    return true;
  }
  else{
    bool flag = true;
    for(auto it = loc_topology_map.begin(); it != loc_topology_map.end(); ++it){
      if(it->first == loc_idx){
        flag = false;
      }
    }
    return flag;
  }

}

//-----------------
// open_ascii_input
//-----------------
void open_ascii_input(void){
    if (!INPUT_FILE) 
    {
      std::cerr << "No input file!" << std::endl;
    }  
    infile = new ifstream(INPUT_FILE);
    if (! infile->is_open()) 
    {
      std::cerr << "Unable to open file \"" << INPUT_FILE << "\" for reading."
      << std::endl;
      exit(-2);
    }
  }

//-----------------
// open_hddm_output
//-----------------
void open_hddm_output(void){
  // Determine output filename from input filename
  OUTPUT_FILE = INPUT_FILE;
  size_t pos = OUTPUT_FILE.find_last_of(".");
  if (pos != string::npos) OUTPUT_FILE.erase(pos);
  OUTPUT_FILE += ".hddm";
  
  // Open output file
  outfile = new ofstream(OUTPUT_FILE.c_str());
  if (! outfile->is_open()) {
      std::cerr << "Unable to open output file \"" << OUTPUT_FILE
                << "\" for writing." << std::endl;
      exit(-3);
   }
  outstream = new hddm_s::ostream(*outfile);

  std::cerr << "Opened output file \"" << OUTPUT_FILE 
            << "\" for writing." << std::endl;
}

//-----------------
// close_hddm_output
//-----------------
void close_files(void){
  // Close input file
  delete infile;

  // Close output file
  delete outstream;
  delete outfile;

  std::cout << "Closed HDDM output file." << std::endl;
}


//-----------------
// clean arrays
//-----------------
void clean_arr(){
  //clear data
  arr_Particletype.clear();
  arr_id.clear();
  arr_parentid.clear();
  arr_E.clear();
  arr_px.clear(); 
  arr_py.clear();
  arr_pz.clear();
  arr_x.clear();
  arr_y.clear();
  arr_z.clear();
  arr_vx.clear();
  arr_vy.clear();
  arr_vz.clear();
  arr_vt.clear();
  arr_t.clear();
}

//-----------------
// open_ascii_input
//-----------------
void parse_kinematic_values_from_ASCII(void){ 

  // parse the basic parameters
  *infile >> nParticles >> beamEnergy >> mechFlag;
  if (beamEnergy == 0 && nParticles == 0){
     std::cerr << "Warning! beamEnergy and nParticles are both zero." << std::endl;
  }
  //read in the entire event
  for(int i = 0; i < nParticles; i++){
      char loc_typestr[256];
      int loc_id, loc_parentid;
      float loc_E, loc_px, loc_py, loc_pz, loc_x, loc_y, loc_z, loc_vx, loc_vy, loc_vz;
      *infile >> loc_id >> loc_parentid >> loc_typestr >> loc_E >> loc_px >> loc_py >> loc_pz;
      *infile >> loc_x >> loc_y >> loc_z >> loc_vx >> loc_vy >> loc_vz;
      //particle type
      int loc_type = Str2GeantParticleID(loc_typestr);
      if(loc_type < 0) loc_type = atoi(loc_typestr);
      //flytime
      float loc_distance = sqrt((loc_x-loc_vx)*(loc_x-loc_vx) + (loc_y-loc_vy)*(loc_y-loc_vy) + (loc_z-loc_vz)*(loc_z-loc_vz));
      float loc_flytime = loc_distance*loc_E/(29.9792458*sqrt(loc_px*loc_px + loc_py*loc_py + loc_pz*loc_pz));
      //store values in the vector
      arr_Particletype.push_back(loc_type);
      arr_id.push_back(loc_id-1);  //start from 1 (in ASCII, mc_gen let it start from 2)
      arr_parentid.push_back(loc_parentid-1); //start from 0 (for the particles from the primary vertex) (in ASCII, mc_gen let it start from 1)
      arr_E.push_back(loc_E);
      arr_px.push_back(loc_px);
      arr_py.push_back(loc_py);
      arr_pz.push_back(loc_pz);
      arr_x.push_back(loc_x);
      arr_y.push_back(loc_y);
      arr_z.push_back(loc_z);
      arr_vx.push_back(loc_vx);
      arr_vy.push_back(loc_vy);
      arr_vz.push_back(loc_vz);
      arr_t.push_back(loc_flytime);
      //cout<<i<<","<<loc_flytime<<endl;
    }
    //cout<<"----------------"<<endl;
}

//-------------------------------
// topology reader
//-------------------------------
void get_topology(void){
    // read information from arr_parentid and arr_id,
    // write the topology in a map<int,set<int>>
    // thus map.size = n_vertices 
    // iterate over the map to fill the event in hddm format

    // clear the memory from last event
        //which is not necessary for current generated events
        // but might be necessary for future situations where multiple reactions are simulated together
    loc_topology_map.clear();

    // loop over id and parent_id to retrieve the topology
    // Example: {1 => { 2, 3, 4 }, 3 => { 5, 6 }, 4 => { 7, 8 }}
    for(int i = 0; i< nParticles; i++){
        loc_topology_map[arr_parentid.at(i)].insert(arr_id.at(i));
    }
}



//-------------------------------
// get_flytime
//-------------------------------
void get_vt(void){
  //calculate and add the flytime recursively
  for(int i = 0; i<nParticles; i++){
    int parent = arr_parentid.at(i);
    float loc_t = 0.0;
    while(parent != 0){
      loc_t+=arr_t.at(parent-1);
      parent = arr_parentid.at(parent-1);
    }
    
    

    arr_vt.push_back(loc_t);
    
  }

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
            case 'r':
               USER_RUNNUMBER = atof(&ptr[1]);
               break;
            default:
              std::cerr << "Unknown option \"" << argv[i] << "\"" << std::endl;
              Usage();
              exit(-1);
         }
      }
      else {
         INPUT_FILE = argv[i];
      }
   }
   

   

}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "       GEN2HDDM [options] file.ascii" << std::endl;
  std::cout << std::endl;
  std::cout << "Convert an ascii file of events generated by mc_gen into HDDM"
            << std::endl;
  std::cout << "for use as input to hdgeant." << std::endl;
  std::cout << std::endl;
  std::cout << " options:" << std::endl;
  std::cout << std::endl;
  std::cout << "  -r#                       "
               "Set the run number (overiding what's in input file)" << std::endl;
  std::cout << "  -h                        "
               "print this usage statement." << std::endl;
  std::cout << std::endl;
}

//-------------------------------
// Str2GeantParticleID
//-------------------------------
int Str2GeantParticleID(char *str)
{
   if (strcmp(str, "unknown") == 0 || strcmp(str, "Unknown") == 0)
      return Unknown;
   if (strcmp(str, "gamma") == 0 || strcmp(str, "Gamma") == 0)
      return Gamma;
   if (strcmp(str, "positron") == 0 || strcmp(str, "Positron") == 0)
      return Positron;
   if (strcmp(str, "electron") == 0 || strcmp(str, "Electron") == 0)
      return Electron;
   if (strcmp(str, "neutrino") == 0 || strcmp(str, "Neutrino") == 0)
      return Neutrino;
   if (strcmp(str, "mu+") == 0 || strcmp(str, "Mu+") == 0)
      return MuonPlus;
   if (strcmp(str, "mu-") == 0 || strcmp(str, "Mu-") == 0)
      return MuonMinus;
   if (strcmp(str, "pi0") == 0 || strcmp(str, "Pi0") == 0)
      return Pi0;
   if (strcmp(str, "pi+") == 0 || strcmp(str, "Pi+") == 0)
      return PiPlus;
   if (strcmp(str, "pi-") == 0 || strcmp(str, "Pi-") == 0)
      return PiMinus;
   if (strcmp(str, "KL") == 0)
      return KLong;
   if (strcmp(str, "K+") == 0)
      return KPlus;
   if (strcmp(str, "K-") == 0)
      return KMinus;
   if (strcmp(str, "neutron") == 0 || strcmp(str, "Neutron") == 0)
      return Neutron;
   if (strcmp(str, "proton") == 0 || strcmp(str, "Proton") == 0)
      return Proton;
   if (strcmp(str, "pbar") == 0 || strcmp(str, "Pbar") == 0)
      return AntiProton;
   if (strcmp(str, "Ks") == 0)
      return KShort;
   if (strcmp(str, "eta") == 0 || strcmp(str, "Eta") == 0)
      return Eta;
   if (strcmp(str, "lambda") == 0 || strcmp(str, "Lambda") == 0)
      return Lambda;
   if (strcmp(str, "sigma+") == 0 || strcmp(str, "Sigma+") == 0)
      return SigmaPlus;
   if (strcmp(str, "sigma0") == 0 || strcmp(str, "Sigma0") == 0)
      return Sigma0;
   if (strcmp(str, "sigma-") == 0 || strcmp(str, "Sigma-") == 0)
      return SigmaMinus;
   if (strcmp(str, "Xi0") == 0)
      return Xi0;
   if (strcmp(str, "Xi-") == 0)
      return XiMinus;
   if (strcmp(str, "omega-") == 0 || strcmp(str, "Omega-") == 0)
      return OmegaMinus;
   if (strcmp(str, "nbar") == 0 || strcmp(str, "Nbar") == 0)
      return AntiNeutron;
   if (strcmp(str, "lambdabar") == 0 || strcmp(str, "Lambdabar") == 0)
      return AntiLambda;
   if (strcmp(str, "sigmabar-") == 0)
      return AntiSigmaMinus;
   if (strcmp(str, "sigmabar0") == 0 || strcmp(str, "Sigmabar0") == 0)
      return AntiSigma0;
   if (strcmp(str, "sigmabar+") == 0 || strcmp(str, "Sigmabar+") == 0)
      return AntiSigmaPlus;
   if (strcmp(str, "Xibar0") == 0)
      return AntiXi0;
   if (strcmp(str, "Xibar+") == 0)
      return AntiXiPlus;
   if (strcmp(str, "omegabar+") == 0 || strcmp(str, "Omegabar+") == 0)
      return AntiOmegaPlus;
   if (strcmp(str, "rho0") == 0 || strcmp(str, "Rho0") == 0)
      return Rho0;
   if (strcmp(str, "rho+") == 0 || strcmp(str, "Rho+") == 0)
      return RhoPlus;
   if (strcmp(str, "rho") == 0 || strcmp(str, "Rho") == 0)
      return RhoMinus;
   if (strcmp(str, "omega") == 0 || strcmp(str, "Omega") == 0)
      return omega;
   if (strcmp(str, "etaprime") == 0 || strcmp(str, "Etaprime") == 0)
      return EtaPrime;
   if (strcmp(str, "phi") == 0 || strcmp(str, "Phi") == 0)
      return phiMeson;
   if (!strcmp(str, "Pb208"))
      return Pb208;
   
   return -1;
}























