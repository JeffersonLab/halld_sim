/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2013-2019  GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Benidikt Zihlmann, Igal Jaegle                            *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#ifndef HDDMOUT_H_
#define HDDMOUT_H_

using namespace std;

#include "HDDM/hddm_s.h"

struct tmpEvt_t {
  int nGen;
  TString str_target;
  TString str_meson;
  TString str_participant;
  TString str_spectator;
  Particle_t t_targ;
  Particle_t t_meso;
  Particle_t t_spec;
  Particle_t t_part;
  double weight;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector q1;
  TLorentzVector q2;
  TLorentzVector q3;
  TLorentzVector q4;
  TLorentzVector q5;
  TLorentzVector q6;
  TLorentzVector q7;
};

class HddmOut {
 private:
  s_iostream_t* ostream;
  //TDatabasePDG* pdg;
  s_PhysicsEvents_t* phyEvt;
  s_Reactions_t* reactions;
  s_Reaction_t* reaction;
  s_Target_t* target;
  s_Beam_t* beam;
  s_Vertices_t* vertices;
  s_HDDM_t* hddmEvt;
  s_Origin_t* origin;
  s_Products_t* products;

  Particle_t targetType;
  Particle_t beamType;
  
 public:
  HddmOut(string filename) {
      cout << "opening HDDM file: " << filename << endl;
    ostream = init_s_HDDM((char*)filename.c_str());
    targetType = Helium;
    beamType = Gamma;
  }
  
  ~HddmOut() {
    close_s_HDDM(ostream);
  }
  
  void init(int runNo) {
    //This sets the run number and event characteristics
    //The HDDM entry has one event, which has one reaction
    hddmEvt = make_s_HDDM();
    hddmEvt->physicsEvents = phyEvt = make_s_PhysicsEvents(1);
    phyEvt->mult = 1;
    phyEvt->in[0].runNo = runNo;
    
    //We define beam and target parameters for the reaction, which
    //remain the same between events
    phyEvt->in[0].reactions = reactions = make_s_Reactions(1);
    reactions->mult = 1;
    reaction = &reactions->in[0];
    reaction->target = target = make_s_Target();
    //target->type = targetType;
    target->properties = make_s_Properties();
    //target->properties->charge = ParticleCharge(targetType);
    //target->properties->mass = ParticleMass(targetType);
    target->momentum = make_s_Momentum();
    //target->momentum->px = 0;
    //target->momentum->py = 0;
    //target->momentum->pz = 0;
    //target->momentum->E  = ParticleMass(targetType);
    reaction->beam = beam = make_s_Beam();
    beam->type = beamType;
    beam->properties = make_s_Properties();
    beam->properties->charge = ParticleCharge(beamType);
    beam->properties->mass = ParticleMass(beamType);
    beam->momentum = make_s_Momentum();

  }
  
  void write(tmpEvt_t evt, int runNum, int eventNum) {
    init(runNum);
    phyEvt->in[0].eventNo = eventNum;
    reaction->vertices = vertices = make_s_Vertices(1);
    vertices->mult = 1;
    vertices->in[0].origin = origin = make_s_Origin();
    vertices->in[0].products = products = make_s_Products(evt.nGen);
    
    origin->t = 0.0;
    origin->vx = 0.0;
    origin->vy = 0.0;
    origin->vz = 0.0;

    target->type = evt.t_targ;
    target->properties->charge = ParticleCharge(evt.t_targ);
    target->properties->mass = ParticleMass(evt.t_targ);
    
    target->momentum->px = evt.target.Px();
    target->momentum->py = evt.target.Py();
    target->momentum->pz = evt.target.Pz();
    target->momentum->E  = evt.target.E();
    
    beam->momentum->px = evt.beam.Px();
    beam->momentum->py = evt.beam.Py();
    beam->momentum->pz = evt.beam.Pz();
    beam->momentum->E  = evt.beam.E();
    //cout <<"beam energy " << evt.beam.E() << " nGen " << evt.nGen << " str_meson " << evt.str_meson << endl;
    products->mult = evt.nGen;
    reaction->weight = evt.weight;
    
    if (evt.nGen == 2) {
      products->in[0].type = evt.t_meso;
      products->in[0].pdgtype = PDGtype(evt.t_meso);
      products->in[0].id = 1;
      products->in[0].parentid = 0;
      products->in[0].mech = 0;
      products->in[0].momentum = make_s_Momentum();
      products->in[0].momentum->px = evt.q1.Px();
      products->in[0].momentum->py = evt.q1.Py();
      products->in[0].momentum->pz = evt.q1.Pz();
      products->in[0].momentum->E = evt.q1.E();
      
      products->in[1].pdgtype = PDGtype(evt.t_targ);
      if (evt.str_target == "Deuteron") {
	//products->in[1].type = Deuteron;
	products->in[1].pdgtype = 1000010020;
      } else if (evt.str_target == "H3" || evt.str_target == "Triton") {
	//products->in[1].type = Triton;
	products->in[1].pdgtype = 1000010030;
      } else if (evt.str_target == "He3" || evt.str_target == "Helium-3") {
	//products->in[1].type = Helium-3;
	products->in[1].pdgtype = 1000020030;
      } else if (evt.str_target == "He4" || evt.str_target == "Helium") {
	//products->in[1].type = Helium;
	products->in[1].pdgtype = 1000020040;
      } else if (evt.str_target == "Be9" || evt.str_target == "Beryllium-9") {
	//products->in[1].type = Be9;
	products->in[1].pdgtype = 1000040090;
      } else if (evt.str_target == "Proton") {
	//products->in[1].type = Proton;
	products->in[1].pdgtype = 2212;
      } else if (evt.str_target == "Neutron") {
	//products->in[1].type = Neutron;
	products->in[1].pdgtype = 2112;
      }
      products->in[1].type = evt.t_targ;
      //cout << PDGtype(evt.t_targ) << endl;
      products->in[1].id = 2;
      products->in[1].parentid = 0;
      products->in[1].mech = 0;
      products->in[1].momentum = make_s_Momentum();
      products->in[1].momentum->px = evt.q2.Px();
      products->in[1].momentum->py = evt.q2.Py();
      products->in[1].momentum->pz = evt.q2.Pz();
      products->in[1].momentum->E = evt.q2.E();
    } else if (evt.nGen == 3) {
    }
    flush_s_HDDM(hddmEvt, ostream);

  }
};

#endif /* HDDMOUT_H_ */
