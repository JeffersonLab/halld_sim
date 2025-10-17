/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#ifndef HDDMOUT_H_
#define HDDMOUT_H_

using namespace std;

#include "HDDM/hddm_s.h"

struct tmpEvt_t {
  int nGen;
  TString rxn;
  double weight;
  Particle_t t_target;
  Particle_t t_spectator;
  Particle_t t_participant;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector q[10];
  Particle_t t_particles[10];
  TLorentzVector recoil;
  int pdg[10];
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
  Particle_t t_target;
  Particle_t t_spectator;
  Particle_t t_participant;
  Particle_t targetType;
  Particle_t beamType;
  
 public:
  HddmOut(string filename) {
      cout << "opening HDDM file: " << filename << endl;
    ostream = init_s_HDDM((char*)filename.c_str());
    targetType = Proton;
    beamType = Gamma;
  }

  Particle_t TYPE;

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
    target->type = targetType;
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

    beam->momentum->px = evt.beam.Px();
    beam->momentum->py = evt.beam.Py();
    beam->momentum->pz = evt.beam.Pz();
    beam->momentum->E  = evt.beam.E();

    target->type = evt.t_target;
    target->properties->charge = ParticleCharge(evt.t_target);
    target->properties->mass = ParticleMass(evt.t_target);
    
    target->momentum->px = evt.target.Px();
    target->momentum->py = evt.target.Py();
    target->momentum->pz = evt.target.Pz();
    target->momentum->E  = evt.target.E();
    
    beam->momentum->px = evt.beam.Px();
    beam->momentum->py = evt.beam.Py();
    beam->momentum->pz = evt.beam.Pz();
    beam->momentum->E  = evt.beam.E();
    
    products->mult = evt.nGen;
    reaction->weight = evt.weight;

    //cout <<"nGen " << evt.nGen << endl;
    
    for (int i = 0; i < evt.nGen; i ++) {
      products->in[i].type = evt.t_particles[i];
      products->in[i].pdgtype = PDGtype(evt.t_particles[i]);
      products->in[i].id = i + 1;
      products->in[i].parentid = 0;
      products->in[i].mech = 0;
      products->in[i].momentum = make_s_Momentum();
      products->in[i].momentum->px = evt.q[i].Px();
      products->in[i].momentum->py = evt.q[i].Py();
      products->in[i].momentum->pz = evt.q[i].Pz();
      products->in[i].momentum->E = evt.q[i].E();
      //products->in[i].momentum->t = evt.q[i].t();
    }
    
    flush_s_HDDM(hddmEvt, ostream);

  }
};

#endif /* HDDMOUT_H_ */
