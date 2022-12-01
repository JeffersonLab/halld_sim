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

#include "UTILITIES/BeamProperties.h"
#ifdef HAVE_EVTGEN
#include "EVTGEN_MODELS/RegisterGlueXModels.h"

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
#endif //EVTGEN_EXTERNAL
#endif //HAVE_EVTGEN

typedef struct {
	bool decayed = false; // not used for now
	int parent_id;
	vector<Particle_t> ids;
	vector<TLorentzVector> p4vs;
} secondary_decay_t;

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

//IA 
TString m_str_Nucleus = "";
TString m_str_Participant = "";
TString m_str_Spectator = "";
TH1F * m_h_PFermi;
double m_mass_nuclei = 0;
double m_ParticipantMass = 0;
double m_SpectatorMass = 0;

// Masses
const double m_p=0.93827; // GeV
const double m_p_sq=m_p*m_p;
double m_eta=0.54775; // GeV
double m_eta_sq=m_eta*m_eta;
// Width
double width=0.;
// Coupling constants 
double g_rho_eta_gamma=0.81;
double g_omega_eta_gamma=0.29;
double g_eta_gamma_gamma=0.0429;
double g_phi_eta_gamma=0.38;

int Nevents=10000;
int runNo=10000;
bool debug=false;

// Diagnostic histograms
TH1D *thrown_t;
TH1D *thrown_dalitzZ;
TH1D *thrown_Egamma;
TH2D *thrown_dalitzXY;  
TH2D *thrown_theta_vs_p;
TH2D *thrown_theta_vs_p_eta;
TH1D *cobrems_vs_E;
TH1D *thrown_FermiP;
TH1D *thrown_f;

char input_file_name[250]="eta548.in";
char output_file_name[250]="eta_gen.hddm";

// Non-default option to generate uniform t-distribution from tmin to tmax
/// (fixed to cross section at tflat_min)
bool gen_uniform_t=false;
float tflat_min=100.; // Physical values are negative. Cross section at larger |t| is equal to cross section at this number (for fixed E_gamma)
float tflat_max=100.; // Physical values are negative

void Usage(void){
  printf("genEtaRegge: generator for eta production based on Regge trajectory formalism.\n");
  printf(" Usage:  genEtaRegge <options>\n");
  printf("   Options:  -N<number of events> (number of events to generate)\n");
  printf("             -O<output.hddm>   (default: eta_gen.hddm)\n");
  printf("             -I<input.in>      (default: eta548.in)\n");
  printf("             -R<run number>    (default: 10000)\n");
  printf("             -h                (Print this message and exit.)\n");
  printf("Coupling constants, photon beam energy range, and eta decay products are\n");
  printf("specified in the <input.in> file.\n");

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
      case 'I':
	sscanf(&ptr[2],"%s",input_file_name);
	break;
      case 'O':
	sscanf(&ptr[2],"%s",output_file_name);
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

// Cross section dsigma/dt derived from Laget(2005)
double CrossSection(double s,double t,double p_gamma,double p_eta,double theta){
  // Coupling constants 
  double c_rho_p_p=0.92/137.;
  double c_omega_p_p=6.44/137.;
  double c_gamma_p_p=1./(137.*137.);
  double c_phi_p_p=0.72/137.;

  double kappa_rho=6.1;
  double kappa_gamma=1.79;
  //double kappa_phi=0.;

  // Angular quantities
  double sintheta=sin(theta);
  double costheta=cos(theta);

  // Kinematic quantities for exchange particle
  double q0=p_gamma-sqrt(p_eta*p_eta+m_eta_sq);
  double q3=p_gamma-p_eta*costheta;
  double q0_minus_q3=q0-q3;
  double q0_minus_q3_sq=q0_minus_q3*q0_minus_q3;
  double q0_sq=q0*q0;
  double q1sq_plus_q2sq=p_eta*p_eta*sintheta*sintheta;
  // Kinematic quantities for target
  double pt3=-p_gamma;
  double pt0=sqrt(m_p_sq+pt3*pt3);
  double pt0_minus_pt3=pt0-pt3;
  double pt0_plus_pt3=pt0+pt3;
  // Kinematic quantities for beam
  double p_gamma_sq=p_gamma*p_gamma;
  // other kinematic quantities
  double pt_dot_q=pt0*q0-pt3*q3;

 // Mass scale for Regge propagators
  double s0=1.0;

  // Regge trajectory for rho
  double a_rho=0.55+0.8*t;
  double a_rho_prime=0.8;
  double regge_rho=pow(s/s0,a_rho-1.)*M_PI*a_rho_prime/(sin(M_PI*a_rho)*TMath::Gamma(a_rho));

  // Regge trajectory for omega
  double a_omega=0.44+0.9*t;
  double a_omega_prime=0.9;
  double regge_omega=pow(s/s0,a_omega-1.)*M_PI*a_omega_prime/(sin(M_PI*a_omega)*TMath::Gamma(a_omega));

  // Regge trajectory for phi(1020)
  double a_phi=0.23+0.7*t;
  double a_phi_prime=0.7;
  double regge_phi=pow(s/s0,a_phi-1.)*M_PI*a_phi_prime/(sin(M_PI*a_phi)*TMath::Gamma(a_phi));

  //  printf("> 36:  %f\n",p_gamma*p_eta*sqrt(mass_factor));

  // amplitude factors for terms involving 0, 1 and 2 powers of kappa
  double amp_factor_kappa0
    =8.*p_gamma_sq*(q1sq_plus_q2sq*(s+q0_minus_q3*pt0_minus_pt3)
		    +q0_minus_q3_sq*pt_dot_q);
  double amp_factor_kappa1=32.*p_gamma_sq*m_p*(q0_minus_q3_sq*t
					       +2.*q0_sq*q1sq_plus_q2sq);
  double amp_factor_kappa2
    =32.*p_gamma_sq*(q1sq_plus_q2sq*(2.*q0_sq*(pt_dot_q-2.*m_p_sq)
				     -t*pt0_plus_pt3*pt0_plus_pt3
				     +4.*pt_dot_q*q0*pt0_plus_pt3)
		     +q0_minus_q3_sq*(t*(pt_dot_q-2.*m_p_sq)
				      +2.*pt_dot_q*pt_dot_q)
		     );

  // rho amplitude 
  double M_rho_sq=16.*M_PI*M_PI*(g_rho_eta_gamma*g_rho_eta_gamma/m_eta_sq)*c_rho_p_p*regge_rho*regge_rho
    *(amp_factor_kappa0-kappa_rho/(4.*m_p)*amp_factor_kappa1
      +kappa_rho*kappa_rho/(16.*m_p_sq)*amp_factor_kappa2);
  double M_rho=-sqrt(M_rho_sq);

  // omega amplitude 
  double M_omega_sq=16.*M_PI*M_PI*(g_omega_eta_gamma*g_omega_eta_gamma/m_eta_sq)*c_omega_p_p*regge_omega*regge_omega
    *amp_factor_kappa0;
  double M_omega=-sqrt(M_omega_sq);

  // phi amplitude 
  double M_phi_sq=16.*M_PI*M_PI*(g_phi_eta_gamma*g_phi_eta_gamma/m_eta_sq)*c_phi_p_p*regge_phi*regge_phi
    *amp_factor_kappa0;
  double M_phi=+sqrt(M_phi_sq);

  // Primakoff amplitude 
  double M_primakoff_sq=16.*M_PI*M_PI*(g_eta_gamma_gamma*g_eta_gamma_gamma/m_eta_sq)*c_gamma_p_p/(t*t)
    *(amp_factor_kappa0-kappa_gamma/(4.*m_p)*amp_factor_kappa1+kappa_gamma*kappa_gamma/(16.*m_p_sq)*amp_factor_kappa2);
  double M_primakoff=sqrt(M_primakoff_sq);

    
  // M_primakoff=0.;
  //M_primakoff_sq=0.;
  
  //M_omega=0.;
  // M_omega_sq=0.;
  
  //M_rho=0.;
  // M_rho_sq=0.;
 
  //M_phi_sq=0.;
  //M_phi=0.;
  
  double pi_a_omega=M_PI*a_omega;
  double pi_a_rho=M_PI*a_rho;
  double pi_a_phi=M_PI*a_phi;
  double M_sq =M_omega_sq+M_rho_sq+M_primakoff_sq+M_phi_sq
    +2.*M_omega*M_phi*cos(pi_a_omega-pi_a_phi)
    +2.*M_omega*M_rho*cos(pi_a_omega-pi_a_rho)
    +2.*M_omega*M_primakoff*cos(pi_a_omega)
    +2.*M_rho*M_primakoff*cos(pi_a_rho)
    +2.*M_rho*M_phi*cos(pi_a_rho-pi_a_phi)
    +2.*M_phi*M_primakoff*cos(pi_a_phi)
    ;
  
  double hbarc_sq=389.; // Convert to micro-barns
  double dsigma_dt=hbarc_sq*M_sq/(4.*64.*M_PI*s*p_gamma_sq);
  // the extra factor for is for 2 photon spins x 2 proton spins

  return(dsigma_dt);
				 

}

// Put particle data into hddm format and output to file
void WriteEvent(unsigned int eventNumber,TLorentzVector &beam, float vert[3],
		vector<Particle_t> &particle_types,
		vector<TLorentzVector> &particle_vectors, 
		vector<bool> &particle_decayed,
		vector< secondary_decay_t > &secondary_vertices,
		s_iostream_t *file){  
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
   if (m_str_Participant == "Proton") 
     ta->type = Proton;
   else if (m_str_Participant == "Neutron")
     ta->type = Neutron;
   else if (m_str_Participant == "")
     ta->type = Proton;
   //ta->type = Proton;
   ta->properties = make_s_Properties();
   ta->properties->charge = ParticleCharge(ta->type);
   ta->properties->mass = ParticleMass(ta->type);
   ta->momentum = make_s_Momentum();
   //ta->momentum->px = 0.;
   //ta->momentum->py = 0.;
   //ta->momentum->pz = 0.;
   //ta->momentum->E  = ParticleMass(ta->type);
   if (m_str_Participant == "" || (m_str_Nucleus == "" && m_str_Participant !=0)) {
     ta->momentum->px = 0.;
     ta->momentum->py = 0.;
     ta->momentum->pz = 0.;
     ta->momentum->E  = ParticleMass(ta->type);
   } else if (m_str_Nucleus != "") {
     ta->momentum->px = target.Px();
     ta->momentum->py = target.Py();
     ta->momentum->pz = target.Pz();
     ta->momentum->E  = target.E();
   }
   // Primary vertex 
   int num_vertices = 1 + secondary_vertices.size();
   rs->in[0].vertices = vs = make_s_Vertices(num_vertices);
   vs->mult = num_vertices;
   vs->in[0].origin = origin = make_s_Origin();
   vs->in[0].products = ps = make_s_Products(particle_vectors.size());
   ps->mult = 0;
   origin->t = 0.0;
   origin->vx = vert[0];
   origin->vy = vert[1];
   origin->vz = vert[2];
   // Final state particles
   int part_ind = 1;
   for (unsigned int i=0;i<particle_vectors.size();i++,ps->mult++){
     Particle_t my_particle=particle_types[i];
     if(particle_decayed[i])
	     ps->in[ps->mult].type = Unknown;  // zero out particle type info so that hdgeant won't decay the particle.  maybe there is a better way?
	 else
	     ps->in[ps->mult].type = my_particle;
     ps->in[ps->mult].pdgtype = PDGtype(my_particle);
     ps->in[ps->mult].id = part_ind++; /* unique value for this particle within the event */
     ps->in[ps->mult].parentid = 0;  /* All internally generated particles have no parent */
     ps->in[ps->mult].mech = 0; // ???     
     ps->in[ps->mult].momentum = make_s_Momentum();
     ps->in[ps->mult].momentum->px = particle_vectors[i].Px();
     ps->in[ps->mult].momentum->py = particle_vectors[i].Py();
     ps->in[ps->mult].momentum->pz = particle_vectors[i].Pz();
     ps->in[ps->mult].momentum->E  = particle_vectors[i].E();
   }
   // write out any secondary vertices (like pi0 decays)
   //		vector< secondary_decay_t > &secondary_vertices,
   if(secondary_vertices.size() > 0) {
   		int vertex_ind = 1;
   		for(auto& the_vertex : secondary_vertices) {
		   // assume that all particles are generated at the same vertex
		   vs->in[vertex_ind].origin = origin = make_s_Origin();
		   vs->in[vertex_ind].products = ps = make_s_Products(the_vertex.ids.size());
		   ps->mult = 0;
		   origin->t = 0.0;
		   origin->vx = vert[0];
		   origin->vy = vert[1];
		   origin->vz = vert[2];
	   
		   // add in the particles associated with this vertex
		   for (unsigned int i=0; i<the_vertex.ids.size(); i++,ps->mult++){
			 Particle_t my_particle = the_vertex.ids[i];
			 ps->in[ps->mult].decayVertex = vertex_ind;
			 ps->in[ps->mult].type = my_particle;
			 ps->in[ps->mult].pdgtype = PDGtype(my_particle);
			 ps->in[ps->mult].id = part_ind++; /* unique value for this particle within the event */
			 ps->in[ps->mult].parentid = the_vertex.parent_id;  /* All internally generated particles have no parent */
			 ps->in[ps->mult].mech = 0; // ???     
			 ps->in[ps->mult].momentum = make_s_Momentum();
			 ps->in[ps->mult].momentum->px = the_vertex.p4vs[i].Px();
			 ps->in[ps->mult].momentum->py = the_vertex.p4vs[i].Py();
			 ps->in[ps->mult].momentum->pz = the_vertex.p4vs[i].Pz();
			 ps->in[ps->mult].momentum->E  = the_vertex.p4vs[i].E();
		   }
		   
		   vertex_ind++;  // get ready for next iteration
	   }
	}
   flush_s_HDDM(thisOutputEvent,file);
}

// Create some diagnostic histograms
void CreateHistograms(string beamConfigFile,int num_decay_particles){

  if(gen_uniform_t) thrown_t=new TH1D("thrown_t","Thrown -t distribution",1000,0.,tflat_max);
  else              thrown_t=new TH1D("thrown_t","Thrown -t distribution",1000,0.,3);
  thrown_t->SetXTitle("-t [GeV^{2}]");
  thrown_Egamma=new TH1D("thrown_Egamma","Thrown E_{#gamma} distribution",
			       1000,0,12.);
  thrown_Egamma->SetTitle("E_{#gamma} [GeV]");
  
  thrown_theta_vs_p=new TH2D("thrown_theta_vs_p","Proton #theta_{LAB} vs. p",
			       200,0,2.,180,0.,90.);
  thrown_theta_vs_p->SetXTitle("p [GeV/c]");
  thrown_theta_vs_p->SetYTitle("#theta [degrees]");
  
  thrown_theta_vs_p_eta=new TH2D("thrown_theta_vs_p_eta","#eta #theta_{LAB} vs. p",
			       120,0,12.,180,0.,180.);
  thrown_theta_vs_p_eta->SetXTitle("p [GeV/c]");
  thrown_theta_vs_p_eta->SetYTitle("#theta [degrees]");
  
  if(num_decay_particles==3) {
      thrown_dalitzZ=new TH1D("thrown_dalitzZ","thrown dalitz Z",110,-0.05,1.05);
      thrown_dalitzXY=new TH2D("thrown_dalitzXY","Dalitz distribution Y vs X",100,-1.,1.,100,-1.,1);
  }
  
  BeamProperties beamProp(beamConfigFile);
  cobrems_vs_E = (TH1D*)beamProp.GetFlux();
}


// Create a graph of the cross section dsigma/dt as a function of -t
void GraphCrossSection(double &xsec_max){
  // beam energy in lab
  double Egamma=cobrems_vs_E->GetBinLowEdge(1); // get from CobremsGenerator histogram

  // CM energy
  double s=m_p*(m_p+2.*Egamma);
  double Ecm=sqrt(s);

  // Momenta of incoming photon and outgoing eta and proton in cm frame
  double p_gamma=(s-m_p_sq)/(2.*Ecm);
  double E_eta=(s+m_eta_sq-m_p_sq)/(2.*Ecm);
  double p_eta=sqrt(E_eta*E_eta-m_eta_sq);
  
  // Momentum transfer t
  double p_diff=p_gamma-p_eta;
  double t0=m_eta_sq*m_eta_sq/(4.*s)-p_diff*p_diff;
  
  double sum=0.;
  double t_old=t0;
  double t_array[10000];
  double xsec_array[10000];
  xsec_max=0.;
  for (unsigned int k=0;k<10000;k++){
    double theta_cm=M_PI*double(k)/10000.;
    double sin_theta_over_2=sin(0.5*theta_cm);
    double t=t0-4.*p_gamma*p_eta*sin_theta_over_2*sin_theta_over_2;
    double xsec=CrossSection(s,t,p_gamma,p_eta,theta_cm);
    if (xsec>xsec_max) xsec_max=xsec;
    
    t_array[k]=-t;
    xsec_array[k]=xsec;

    sum-=xsec*(t-t_old);
    t_old=t;
  }
  TGraph *Gxsec=new TGraph(10000,t_array,xsec_array);
  TString xsec_title = "#eta Cross Section at E_{#gamma}="+to_string(Egamma)+" GeV;-t [GeV^{2}];d#sigma/dt [#mub/GeV^{2}]";
  Gxsec->SetTitle(xsec_title);
  Gxsec->Write("Cross section");
 
  cout << "Total cross section at " << Egamma << " GeV = "<< sum 
       << " micro-barns"<<endl;
}


//-----------
// main
//-----------
int main(int narg, char *argv[])
{  
  ParseCommandLineArguments(narg, argv);


  // open ROOT file
  string rootfilename="eta_gen.root";
  TFile *rootfile=new TFile(rootfilename.c_str(),"RECREATE",
			    "Produced by genEta");

  // open HDDM file
  s_iostream_t *file = init_s_HDDM(output_file_name);
 
 
  // Initialize random number generator
  TRandom3 *myrand=new TRandom3(0);// If seed is 0, the seed is automatically computed via a TUUID object, according to TRandom3 documentation

  // Fixed target
  TLorentzVector target(0.,0.,0.,m_p);

  //----------------------------------------------------------------------------
  // Get production (Egamma range) and decay parameters from input file
  //----------------------------------------------------------------------------

  // Start reading the input file 
  ifstream infile(input_file_name);
  if (!infile.is_open()){
    cerr << "Input file missing! Exiting..." <<endl;
    exit(-1);
  } 
  
  // IA, get generator config file
  MyReadConfig * ReadFile = new MyReadConfig();
  ReadFile->ReadConfigFile(input_file_name);
  m_str_Nucleus = ReadFile->GetConfigName("Nucleus");
  m_str_Participant = ReadFile->GetConfigName("Participant");
  m_str_Spectator = ReadFile->GetConfigName("Spectator");
  TString m_str_Fermi_file = ReadFile->GetConfigName("FermiMotionFile");
  if (m_str_Nucleus != "") { 
    if (m_str_Nucleus == "D2") { 
      m_mass_nuclei = 1.875613;
      if(m_str_Participant == "Proton") {
	m_ParticipantMass = 0.93827;
	m_SpectatorMass = 0.93956;
      } else if (m_str_Participant == "Neutron") {
	m_ParticipantMass = 0.93956;
	m_SpectatorMass = 0.93827;
      }
    } else if (m_str_Nucleus == "He4") {
      m_mass_nuclei = 3.727379;
      if(m_str_Participant == "Proton") {
	m_ParticipantMass = 0.93827;
	m_SpectatorMass = 2.808921;//003 001 3H
      } else if (m_str_Participant == "Neutron") {
	m_ParticipantMass = 0.93956;
	m_SpectatorMass = 2.808391;//003 002 3He
      }
    } else if (m_str_Nucleus == "C12") {
      m_mass_nuclei = 11.174862;
      if(m_str_Participant == "Proton") {
	m_ParticipantMass = 0.93827;
	m_SpectatorMass = 10.252547;//011 005 11B
      } else if (m_str_Participant == "Neutron") {
	m_ParticipantMass = 0.93956;
	m_SpectatorMass = 10.254018;//011 006 11C
      }
    }
    m_p = m_ParticipantMass;
    m_p_sq = m_p * m_p;
  } else if (m_str_Nucleus == "" && m_str_Participant != "") {
    if(m_str_Participant == "Proton") {
      m_ParticipantMass = 0.93827;
    } else if (m_str_Participant == "Neutron") {
      m_ParticipantMass = 0.93956;
    }
    m_p = m_ParticipantMass;
    m_p_sq = m_p * m_p;
  }
  if (m_str_Fermi_file != "") {
    cout <<"Target is made of " << m_str_Nucleus << " with the participant " << m_str_Participant << " and spectator " << m_str_Spectator << endl;
    cout <<"Nucleon Fermi motion is located in " << m_str_Fermi_file << endl;
    m_h_PFermi = new TH1F("PFermi", "", 1000, 0.0, 1.0);
    ifstream in;
    in.open(m_str_Fermi_file);
    int i = 0;
    while (in.good()) {
      double pf = 0, val = 0;
      in >> pf >> val;
      if (val > 0) {
	m_h_PFermi->SetBinContent(i + 1, val);
	i ++;
      }
    }
    in.close();
  }

  // Get beam properties configuration file
  string comment_line;
  getline(infile,comment_line);
  string beamConfigFile;
  infile >> beamConfigFile;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "Photon beam configuration file " << beamConfigFile.data() << endl;

  // Get decaying particle mass and width
  string comment_line2;
  getline(infile,comment_line);
  double m_eta_R=0.;
  infile >> m_eta_R;
  infile >> width;
  infile.ignore(); // ignore the '\n' at the end of this line

  m_eta=m_eta_R;
  m_eta_sq=m_eta*m_eta;

  cout << "Mass, width of decaying particle [GeV] = "<< m_eta <<"," << width << endl;

  // Get coupling constants for photon vertex
  getline(infile,comment_line);
  infile >> g_eta_gamma_gamma;
  infile >> g_rho_eta_gamma;
  infile >> g_omega_eta_gamma;
  infile >> g_phi_eta_gamma;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "Coupling constants:" <<endl;
  cout << " g_eta_gamma_gamma = " << g_eta_gamma_gamma <<endl;
  cout << " g_rho_eta_gamma = " << g_rho_eta_gamma <<endl;
  cout << " g_omega_eta_gamma = " << g_omega_eta_gamma << endl;
  cout << " g_phi_eta_gamma = " << g_phi_eta_gamma <<endl;
  
  // Get number of decay particles
  int num_decay_particles=0;
  getline(infile,comment_line);
  infile >> num_decay_particles;
  infile.ignore(); // ignore the '\n' at the end of this line

  cout << "number of decay particles = " << num_decay_particles << endl;

  

  bool use_evtgen = false;
#ifdef HAVE_EVTGEN
  // check to see if we should use EvtGen
  EvtParticle* parent(0);
  EvtGen *myGenerator = nullptr;
  EvtId EtaId;  // read this in from file below - this is actually a class, which gets filled with reasonable defaults
  
  // if we don't explicitly define the decay particles in the config file, 
  // then use EvtGen to generate the decays
  if(num_decay_particles <= 0) { 
  	num_decay_particles = 0;
  	use_evtgen = true;

    cout << "Using EvtGen to decay particles ..." << endl;

	// get the produced particle type
	getline(infile,comment_line);
	string particle_type;
	getline(infile,particle_type);
	trim(particle_type);
	cout << "Generating particle: " << particle_type << endl;
  	
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

	// Now that we have initialized EvtGen, we can access things like the particle DB
	EtaId = EvtPDL::getId(std::string(particle_type));
  	
  }
#endif // HAVE_EVTGEN
  
  // Set up vectors of particle ids
  vector<Particle_t>particle_types;
  //double *decay_masses =new double[num_decay_particles];
  vector<double> decay_masses(num_decay_particles);
  double *res_decay_masses=NULL;
  vector<Particle_t>res_particle_types;

  // GEANT ids of decay particles
  getline(infile,comment_line);
  cout << comment_line << endl;
  int reson_index=-1;
  if(!use_evtgen) {
	  cout << "Particle id's of decay particles =";
	  for (int k=0;k<num_decay_particles;k++){
		int ipart;
		infile >> ipart;
		cout << " " << ipart; 
		particle_types.push_back((Particle_t)ipart);
		if (ipart>0){
		  decay_masses[k]=ParticleMass((Particle_t)ipart);
		}
		else {
		  reson_index=k;
		}
	  } 
	  cout << endl;
  }
  unsigned int num_res_decay_particles=0; 
  double reson_mass=0.,reson_width=0.;
  int reson_L=0;
  if (reson_index>=0){
    infile.ignore(); // ignore the '\n' at the end of this line 
    getline(infile,comment_line);
    cout << comment_line << endl;
    infile >> reson_mass;
    decay_masses[reson_index]=reson_mass;
    cout << "Resonance mass = " << reson_mass << " [GeV]" << endl;   
    infile >> reson_width;
    cout << "Resonance width = " << reson_width << " [GeV]" << endl; 
    infile >> reson_L;
    cout << "Resonance orbital angular momentum L = " << reson_L << endl; 
    infile >> num_res_decay_particles;
    if (num_res_decay_particles>1) res_decay_masses=new double[num_res_decay_particles]; 
    else{
      cout << "Invalid number of decay particles! " << endl;
      exit(0);
    }
    cout << " Decay particles: ";
    for (unsigned int i=0;i<num_res_decay_particles;i++){
      int ipart;
      infile >> ipart;
      cout << " " << ipart; 
      res_particle_types.push_back((Particle_t)ipart);
      res_decay_masses[i]=ParticleMass((Particle_t)ipart);
    }
  }

  // Search for lines in input file starting with "tflat_min" or "tflat_max", if found we reset globals
  while( !infile.eof() ) {
    string line   = "";
    string tflat_string = "";
    getline(infile,line);
    if(line.length() < 11) continue;
    // Yes this code is ugly, but works. I miss python.
    if(line.substr(0,9) == "tflat_min") {
        string str_tval="";
        for(size_t loc_i=9; loc_i<line.length(); loc_i++) {
            string this_char; this_char += line[loc_i];
            if(isdigit(line[loc_i]) || this_char=="." ) {
                str_tval+=line[loc_i];
            }
            if(str_tval.length()>0 && this_char==" ") break;
        }
        tflat_min = -1*fabs(atof(str_tval.c_str()));
    }
    // Yes this code is ugly, but works. I miss python.
    if(line.substr(0,9) == "tflat_max") {
        string str_tval="";
        for(size_t loc_i=9; loc_i<line.length(); loc_i++) {
            string this_char; this_char += line[loc_i];
            if(isdigit(line[loc_i]) || this_char=="." ) {
                str_tval+=line[loc_i];
            }
            if(str_tval.length()>0 && this_char==" ") break;
        }
        tflat_max = -1*fabs(atof(str_tval.c_str()));
    }
  }
  if(tflat_min<0.&&tflat_max<0.&&tflat_max<tflat_min) gen_uniform_t = true;
  if(gen_uniform_t) cout << "GENERATING DATA WITH UNIFORM T-DIST FROM " << tflat_min << " TO " << tflat_max << endl;
  

  infile.close();
  
  // Create some diagonistic histographs
  CreateHistograms(beamConfigFile,num_decay_particles);

  // Make a TGraph of the cross section at a fixed beam energy
  double xsec_max=0.;
  GraphCrossSection(xsec_max);

  //----------------------------------------------------------------------------
  // Event generation loop
  //----------------------------------------------------------------------------
  for (int i=1;i<=Nevents;i++){
    double Egamma=0.;
    double xsec=0.,xsec_test=0.;

    // Polar angle in center of mass frame
    double theta_cm=0.;

    // Eta momentum in cm
    double p_eta=0.;

    // Transfer 4-momentum;
    double t=0.;

    // vertex position at target
    float vert[4]={0.,0.,0.,0.};

    // IA variables
    double p_Fermi = 0, p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
    double ParticipantEnergy = 0;
    TLorentzVector Ptotal_4Vec(0, 0, 0, 0);
    if (m_str_Nucleus != "") {
      p_Fermi = m_h_PFermi->GetRandom();
      thrown_FermiP->Fill(p_Fermi);
      p_Fermi_x = 0, p_Fermi_y = 0, p_Fermi_z = 0;
      gRandom->Sphere(p_Fermi_x, p_Fermi_y, p_Fermi_z, p_Fermi);
      ParticipantEnergy = m_mass_nuclei - sqrt(pow(m_SpectatorMass, 2) + pow(p_Fermi, 2));
    }

    // use the rejection method to produce eta's based on the cross section
    do{
      // First generate a beam photon using bremsstrahlung spectrum
      Egamma = cobrems_vs_E->GetRandom();

      // CM energy
      double s=m_p*(m_p+2.*Egamma);
      double Ecm=sqrt(s);

      // IA, momenta of incoming photon and outgoing eta and proton in cm frame
      if (m_str_Nucleus != "") {
	Ptotal_4Vec = TLorentzVector(p_Fermi_x, p_Fermi_y, Egamma + p_Fermi_z, Egamma + ParticipantEnergy);
	Ecm = Ptotal_4Vec.M();
	s = pow(Ecm, 2);
      }

      // Momenta of incoming photon and outgoing eta and proton in cm frame
      double p_gamma=(s-m_p_sq)/(2.*Ecm);

      // Generate mass distribution for unstable particle in the final state 
      // with non-negligible width if specified in the input file
      if(!use_evtgen) {
      double mass_check=0.;
      do {
	if (reson_index>-1 && reson_width>0.){
	  if (num_res_decay_particles==2){
	    double BW=0.,BWtest=0.;
	    double m1sq=res_decay_masses[0]*res_decay_masses[0];
	    double m2sq=res_decay_masses[1]*res_decay_masses[1];
	    double m0sq=reson_mass*reson_mass;
	    double m1sq_minus_m2sq=m1sq-m2sq;
	    double q0sq=(m0sq*m0sq-2.*m0sq*(m1sq+m2sq)
			 +m1sq_minus_m2sq*m1sq_minus_m2sq)/(4.*m0sq);
	    double BlattWeisskopf=1.;
	    double d=5.; // Meson radius in GeV^-1: 5 GeV^-1 -> 1 fm
	    double dsq=d*d;
	    double Gamma0sq=reson_width*reson_width;
	    double BWmax=0.;
	    double BWmin=0.;
	    double m_min=res_decay_masses[0]+res_decay_masses[1];
	    //	    double BW_at_m_min=0.;
	    //double BW_at_m_max=0.;
	    if (reson_L==0){
	      BWmax=1./(m0sq*Gamma0sq);
	    }
	    else{
	      BWmax=pow(q0sq,2*reson_L)/(m0sq*Gamma0sq);
	    }
	    double m_max=m_p*(sqrt(1.+2.*Egamma/m_p)-1.);
	    for (int im=0;im<num_decay_particles;im++){
	      if (im==reson_index) continue;
	      m_max-=decay_masses[im];
	    }
	    double m=0.;
	    do {
	      m=m_min+myrand->Uniform(m_max-m_min);
	      double msq=m*m;
	      double qsq=(msq*msq-2.*msq*(m1sq+m2sq)
			  +m1sq_minus_m2sq*m1sq_minus_m2sq)/(4.*msq);
	      double z=dsq*qsq;
	      double z0=dsq*q0sq;
	      double m0sq_minus_msq=m0sq-msq;
	      if (reson_L==0){
		BW=1./(m0sq_minus_msq*m0sq_minus_msq+m0sq*qsq/q0sq*Gamma0sq);
	      }
	      else{
		if (reson_L==1){
		  BlattWeisskopf=(1.+z0)/(1.+z);
		}
		BW=pow(qsq,2*reson_L)
		  /(m0sq_minus_msq*m0sq_minus_msq
		    +m0sq*pow(qsq/q0sq,2*reson_L+1)*Gamma0sq*BlattWeisskopf);
	      }
	      BWtest=BWmin+myrand->Uniform(BWmax-BWmin);
	    } while (BWtest>BW);
	    decay_masses[reson_index]=m;
	  }	
	}

	if (width>0.001){  
	  // Take into account width of resonance, but apply a practical minimum
	  // for the width, overwise we are just wasting cpu cycles...
	  // Use a relativistic Breit-Wigner distribution for the shape.  
	  double m_max_=m_p*(sqrt(1.+2.*Egamma/m_p)-1.);
	  double m_min_=decay_masses[0]; // will add the second mass below
	  double m1sq_=decay_masses[0]*decay_masses[0];
	  double m2sq_=0.;
	  switch(num_decay_particles){
	  case 2:
	    {
	      m_min_+=decay_masses[1];
	      m2sq_=decay_masses[1]*decay_masses[1];
	      break;
	    }
	  case 3:
	    // Define an effective mass in an ad hoc way: we assume that in the 
	    // CM one particle goes in one direction and the two other particles
	    // go in the opposite direction such that p1=-p2-p3.  The effective
	    // mass of the 2-3 system must be something between min=m2+m3 
	    // and max=M-m1, where M is the mass of the resonance.  For
	    // simplicity use the average of these two extremes.
	    {
	      double m2_=0.5*(m_eta_R-decay_masses[0]+decay_masses[1]+decay_masses[2]);
	      m_min_+=decay_masses[1]+decay_masses[2];
	      m2sq_=m2_*m2_;
	      break;
	    }
	  default:
	    break;
	  }
	  double m0sq_=m_eta_R*m_eta_R;
	  double BW_=0.,BWtest_=0.;
	  double Gamma0sq_=width*width;
	  double m1sq_minus_m2sq_=m1sq_-m2sq_;
	  double q0sq_=(m0sq_*m0sq_-2.*m0sq_*(m1sq_+m2sq_)
			+m1sq_minus_m2sq_*m1sq_minus_m2sq_)/(4.*m0sq_);
	  double BWmax_=1./(Gamma0sq_*m0sq_);
	  double BWmin_=0.;
	  double m_=0.;
	  do{
	    m_=m_min_+myrand->Uniform(m_max_-m_min_);
	    double msq_=m_*m_;
	    double qsq_=(msq_*msq_-2.*msq_*(m1sq_+m2sq_)
			 +m1sq_minus_m2sq_*m1sq_minus_m2sq_)/(4.*msq_);
	    double m0sq_minus_msq_=m0sq_-msq_;
	    BW_=1./(m0sq_minus_msq_*m0sq_minus_msq_
		    +m0sq_*m0sq_*Gamma0sq_*qsq_/q0sq_);
	    BWtest_=BWmin_+myrand->Uniform(BWmax_-BWmin_);
	  }
	  while (BWtest_>BW_);
	  m_eta=m_;
	  m_eta_sq=m_*m_;
	}
	// Check that the decay products are consistent with a particle of mass
	// m_eta...
	mass_check=decay_masses[0];
	for (int im=1;im<num_decay_particles;im++){
	  mass_check+=decay_masses[im];
	}
      } while (mass_check>m_eta);
    } else {
    	// if using EvtGen, we don't know what the decay products are a priori
    	// use a non-relativistic BW instead
    	if(width > 0.000001) {
			bool is_good_mass = false;
			do {
				double m_ = myrand->BreitWigner(m_eta_R,width);
				cout << m_ << " " << m_eta_R << " " << width << endl;
				// use a cutoff of +- 5 times the width so we don't populate the far tails too much
				is_good_mass = (m_ > m_eta_R-5.*width) && (m_ < m_eta_R+5.*width);
			} while (!is_good_mass);
		} else {
			m_eta = m_eta_R;
		}
    }
    
    
      double E_eta=(s+m_eta_sq-m_p_sq)/(2.*Ecm);
      p_eta=sqrt(E_eta*E_eta-m_eta_sq);
    
      // Momentum transfer t
      double p_diff=p_gamma-p_eta;
      double t0=m_eta_sq*m_eta_sq/(4.*s)-p_diff*p_diff;
      double sin_theta_over_2=0.;
      t=t0;
      
      // Generate cos(theta) with a uniform distribution and compute the cross 
      // section at this value
      double cos_theta_cm=-1.0+myrand->Uniform(2.);
      theta_cm=acos(cos_theta_cm);
      
      sin_theta_over_2=sin(0.5*theta_cm);
      t=t0-4.*p_gamma*p_eta*sin_theta_over_2*sin_theta_over_2;
      xsec=CrossSection(s,t,p_gamma,p_eta,theta_cm);	  
	  
	  // If generating a sample uniform in t, we need to fix t and re-calculate theta_cm based on it. Others do not depend on t.
	  if(gen_uniform_t&&t<tflat_min&&t>tflat_max) {
		  //Cross section at fixed t value (tflat_min)
		  double t_tmp = tflat_min;
          // if(t0<t_tmp && t_tmp < 0. ) t_tmp=t0-0.00001; //If tflat_min is unphysically small, use (essentially) t0 instead
		  double theta_cm_tmp = 2.*asin(0.5*sqrt( (t0-t_tmp)/(p_gamma*p_eta) ) );
		  xsec=CrossSection(s,t_tmp,p_gamma,p_eta,theta_cm_tmp);
		  // Make t uniform, calculate theta_cm based off of it
		  t=myrand->Uniform(tflat_max,  min( float(t0),tflat_min) ); // If t_min_uniform provided is unphysical, then use physical t_min.
		  theta_cm=2.*asin(0.5*sqrt( (t0-t)/(p_gamma*p_eta) ) );
		  if( std::isnan(theta_cm)==true ) xsec=-1.; // Lazy person's way of skipping unphysical theta_cm. Breaking do/while to accept event will never be satisfied for this case.
	  }

      // Generate a test value for the cross section
      xsec_test=myrand->Uniform(xsec_max);
    }
    while (xsec_test>xsec);
    
    // Generate phi using uniform distribution
    double phi_cm=myrand->Uniform(2.*M_PI);

    // beam 4-vector (ignoring px and py, which are extremely small)
    TLorentzVector beam(0.,0.,Egamma,Egamma);
    thrown_Egamma->Fill(Egamma);

    // IA nucleon/nuclei target
    if (m_str_Nucleus != "") {
      target = TLorentzVector(p_Fermi_x, p_Fermi_y, p_Fermi_z, ParticipantEnergy);
      //cout <<"rewrite target p4 w/ fermi"<<endl;
    } else if (m_str_Nucleus == "" && m_str_Participant != "") {
      target = TLorentzVector(0, 0, 0, m_p);
    }

    // Velocity of the cm frame with respect to the lab frame
    TVector3 v_cm=(1./(Egamma+m_p))*beam.Vect();
    // Four-moementum of the eta in the CM frame
    double pt=p_eta*sin(theta_cm);
    TLorentzVector eta4(pt*cos(phi_cm),pt*sin(phi_cm),p_eta*cos(theta_cm),
			sqrt(p_eta*p_eta+m_eta_sq));

    //Boost the eta 4-momentum into the lab
    //eta4.Boost(v_cm);
    // IA modified boost
    if (m_str_Nucleus != "") { 
      eta4.Boost(Ptotal_4Vec.BoostVector());
    } else if (m_str_Nucleus == "" && m_str_Participant != "") { 
      eta4.Boost(v_cm);
    }


    // Compute the 4-momentum for the recoil proton
    TLorentzVector proton4=beam+target-eta4; 

    //proton4.Print();
    thrown_theta_vs_p->Fill(proton4.P(),180./M_PI*proton4.Theta());
    thrown_theta_vs_p_eta->Fill(eta4.P(),180./M_PI*eta4.Theta());

    // Other diagnostic histograms
    thrown_t->Fill(-t);

    // Gather the particles in the reaction and write out event in hddm format
    vector<TLorentzVector>output_particle_vectors;
    output_particle_vectors.push_back(proton4);
    
    vector<Particle_t>output_particle_types;
    //output_particle_types.push_back(Proton);
    // IA modifications
    if (m_str_Participant == "Proton") 
      output_particle_types.push_back(Proton);
    else if (m_str_Participant == "Neutron") 
      output_particle_types.push_back(Neutron);
    else if (m_str_Participant == "")
      output_particle_types.push_back(Proton);

    vector<bool>output_particle_decays;
    output_particle_decays.push_back(false);
	
    vector<secondary_decay_t>secondary_vertices;
#ifdef HAVE_EVTGEN
    if(use_evtgen) {
        // Set up the parent particle
        EvtVector4R pInit(eta4.E(), eta4.X(), eta4.Y(), eta4.Z());
        parent = EvtParticleFactory::particleFactory(EtaId, pInit);

        if(num_decay_particles == 0) {  
            // just generate the eta and let the decay be done by an external program like decay_evtgen
            // this plays better with MCWrapper, apparently...
            TLorentzVector vec4v( eta4.X(), eta4.Y(), eta4.Z(), eta4.E());

            output_particle_vectors.push_back(vec4v);
            output_particle_types.push_back(PDGtoPType(parent->getPDGId()));
        } else {            
            // Generate the event
            myGenerator->generateDecay(parent);
            
            // Write out resulting particles
            for(unsigned int i=0; i<parent->getNDaug(); i++) {
                TLorentzVector vec4v( parent->getDaug(i)->getP4Lab().get(1),
                                      parent->getDaug(i)->getP4Lab().get(2),
                                      parent->getDaug(i)->getP4Lab().get(3),
                                      parent->getDaug(i)->getP4Lab().get(0)   );
                output_particle_vectors.push_back(vec4v);
                output_particle_types.push_back(PDGtoPType(parent->getDaug(i)->getPDGId()));
						
                // see if any of the particles decay and add info on them
                // should be mostly pi0's, but we should go recursive...
                if(parent->getDaug(i)->getNDaug()>0) {
                    output_particle_decays.push_back(true);
                    
                    secondary_decay_t secondary_vertex;
                    secondary_vertex.parent_id = i;
                    for(unsigned int j=0; j<parent->getDaug(i)->getNDaug(); j++) {
                        TLorentzVector vec4v( parent->getDaug(i)->getDaug(j)->getP4Lab().get(1),
                                              parent->getDaug(i)->getDaug(j)->getP4Lab().get(2),
                                              parent->getDaug(i)->getDaug(j)->getP4Lab().get(3),
                                              parent->getDaug(i)->getDaug(j)->getP4Lab().get(0)   );
                        secondary_vertex.p4vs.push_back(vec4v);
                        secondary_vertex.ids.push_back(PDGtoPType(parent->getDaug(i)->getDaug(j)->getPDGId()));
                    }
                    secondary_vertices.push_back(secondary_vertex);
                } else {
                    output_particle_decays.push_back(false);
                }
            }
	
            parent->deleteTree();
        }
    } else {   // no evtgen
#endif //HAVE_EVTGEN
		// Generate 3-body decay of eta according to phase space
		TGenPhaseSpace phase_space;
		phase_space.SetDecay(eta4,num_decay_particles,decay_masses.data());
		double weight=0.,rand_weight=1.;
		do{
		  weight=phase_space.Generate();
		  rand_weight=myrand->Uniform(1.);
		}
		while (rand_weight>weight);

		// Histograms of Dalitz distribution
		if (num_decay_particles==3){
		  TLorentzVector one=*phase_space.GetDecay(0);  
		  TLorentzVector two=*phase_space.GetDecay(1);
		  TLorentzVector three=*phase_space.GetDecay(2);
		  TLorentzVector one_two=one+two;
		  TLorentzVector one_three=one+three;
 
		  TLorentzVector eta=one_two+three;
		  TVector3 boost=-eta.BoostVector();

		  double eta_mass=eta.M();
		  eta.Boost(boost);
	  
		  one.Boost(boost);
		  two.Boost(boost);
		  three.Boost(boost);

		  double m1=one.M(),m2=two.M(),m3=three.M();
		  double E1=one.E(),E2=two.E(),E3=three.E();
		  double T1=E1-m1; // pi0 for charged channel
		  double T2=E2-m2; // pi+ for charged channel
		  double T3=E3-m3; // pi- for charged channel
		  double Q_eta=eta_mass-m1-m2-m3;
		  double X=sqrt(3.)*(T2-T3)/Q_eta;
		  double Y=3.*T1/Q_eta-1.;
		  thrown_dalitzXY->Fill(X,Y);

		  double z_dalitz=X*X+Y*Y;
		  //printf("z %f\n",z_dalitz);
		  thrown_dalitzZ->Fill(z_dalitz);
		}

		for (int j=0;j<num_decay_particles;j++){
			if (particle_types[j]!=Unknown) {
				output_particle_vectors.push_back(*phase_space.GetDecay(j));
				output_particle_types.push_back(particle_types[j]);
			} else {
				TGenPhaseSpace phase_space2;
				phase_space2.SetDecay(*phase_space.GetDecay(j),num_res_decay_particles,
						  res_decay_masses);
				weight=0.,rand_weight=1.;
				do{
					weight=phase_space2.Generate();
					rand_weight=myrand->Uniform(1.);
				} while (rand_weight>weight);
				for (unsigned int im=0;im<num_res_decay_particles;im++){
					output_particle_types.push_back(res_particle_types[im]);
					output_particle_vectors.push_back(*phase_space2.GetDecay(im));
				}
			}
		}
#ifdef HAVE_EVTGEN
    }
#endif //HAVE_EVTGEN
    
    // Write Event to HDDM file
    WriteEvent(i,beam,vert,output_particle_types,output_particle_vectors,output_particle_decays,secondary_vertices,file);
    
    if (((10*i)%Nevents)==0) cout << 100.*double(i)/double(Nevents) << "\% done" << endl;
  }


  // Write histograms and close root file
  rootfile->Write();
  rootfile->Close();

  // Close HDDM file
  close_s_HDDM(file);
  cout<<endl<<"Closed HDDM file"<<endl;
  cout<<" "<<Nevents<<" event written to "<<output_file_name<<endl;

  // Cleanup
  //delete []decay_masses;
  if (res_decay_masses!=NULL) delete []res_decay_masses;

  return 0;
}
