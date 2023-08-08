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
    target->type = targetType;
    target->properties = make_s_Properties();
    target->properties->charge = ParticleCharge(targetType);
    target->properties->mass = ParticleMass(targetType);
    target->momentum = make_s_Momentum();
    target->momentum->px = 0;
    target->momentum->py = 0;
    target->momentum->pz = 0;
    target->momentum->E  = ParticleMass(targetType);
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
    //cout <<"beam energy " << evt.beam.E() << " nGen " << evt.nGen << " str_meson " << evt.str_meson << endl;
    products->mult = evt.nGen;
    reaction->weight = evt.weight;
    if (evt.nGen == 2) {
      //PRODUCED PHOTON
      if (evt.str_meson == "eta") {
	products->in[0].type = Eta;
	products->in[0].pdgtype = 221;
      } else if (evt.str_meson == "eta'") {
	products->in[0].type = EtaPrime;
	products->in[0].pdgtype = 331;
	//cout <<"eta' " << endl; 
      } else if (evt.str_meson == "pi0") {
	products->in[0].type = Pi0;
	products->in[0].pdgtype = 111;
      }
      products->in[0].id = 1;
      products->in[0].parentid = 0;
      products->in[0].mech = 0;
      products->in[0].momentum = make_s_Momentum();
      products->in[0].momentum->px = evt.q1.Px();
      products->in[0].momentum->py = evt.q1.Py();
      products->in[0].momentum->pz = evt.q1.Pz();
      products->in[0].momentum->E = evt.q1.E();
            
      //PRODUCED Nucleus recoil
      if (evt.str_target == "He4" || evt.str_target == "Helium") {
	products->in[1].type = Helium;
	products->in[1].pdgtype = 1000020040;
      } else if (evt.str_target == "Be9" || evt.str_target == "Beryllium-9") {
	products->in[1].type = Be9;
	products->in[1].pdgtype = 1000040090;
      } else if (evt.str_target == "Proton") {
	products->in[1].type = Proton;
	products->in[1].pdgtype = 2212;
      } else if (evt.str_target == "Neutron") {
	products->in[1].type = Neutron;
	products->in[1].pdgtype = 2112;
      }
      products->in[1].id = 2;
      products->in[1].parentid = 0;
      products->in[1].mech = 0;
      products->in[1].momentum = make_s_Momentum();
      products->in[1].momentum->px = evt.q2.Px();
      products->in[1].momentum->py = evt.q2.Py();
      products->in[1].momentum->pz = evt.q2.Pz();
      products->in[1].momentum->E = evt.q2.E();
    } else if (evt.nGen == 3) {
      //PRODUCED PHOTON
      products->in[0].type = Gamma;
      products->in[0].pdgtype = 22;
      products->in[0].id = 1;
      products->in[0].parentid = 0;
      products->in[0].mech = 0;
      products->in[0].momentum = make_s_Momentum();
      products->in[0].momentum->px = evt.q1.Px();
      products->in[0].momentum->py = evt.q1.Py();
      products->in[0].momentum->pz = evt.q1.Pz();
      products->in[0].momentum->E = evt.q1.E();
      
      //PRODUCED PHOTON
      products->in[1].type = Gamma;
      products->in[1].pdgtype = 22;
      products->in[1].id = 2;
      products->in[1].parentid = 0;
      products->in[1].mech = 0;
      products->in[1].momentum = make_s_Momentum();
      products->in[1].momentum->px = evt.q2.Px();
      products->in[1].momentum->py = evt.q2.Py();
      products->in[1].momentum->pz = evt.q2.Pz();
      products->in[1].momentum->E = evt.q2.E();
      
      //PRODUCED Nucleus recoil
      //products->in[2].type = Helium;
      //products->in[2].pdgtype = 1000020040;
      if (evt.str_target == "He4" || evt.str_target == "Helium") {
	products->in[2].type = Helium;
	products->in[2].pdgtype = 1000020040;
      } else if (evt.str_target == "Be9" || evt.str_target == "Beryllium-9") {
	products->in[2].type = Be9;
	products->in[2].pdgtype = 1000040090;
      } else if (evt.str_target == "Proton") {
	products->in[2].type = Proton;
	products->in[2].pdgtype = 2212;
      } else if (evt.str_target == "Neutron") {
	products->in[2].type = Neutron;
	products->in[2].pdgtype = 2112;
      }
      products->in[2].id = 3;
      products->in[2].parentid = 0;
      products->in[2].mech = 0;
      products->in[2].momentum = make_s_Momentum();
      products->in[2].momentum->px = evt.q3.Px();
      products->in[2].momentum->py = evt.q3.Py();
      products->in[2].momentum->pz = evt.q3.Pz();
      products->in[2].momentum->E = evt.q3.E();
    } else if (evt.nGen == 7) {
      //PRODUCED PHOTON
      products->in[0].type = Gamma;
      products->in[0].pdgtype = 22;
      products->in[0].id = 1;
      products->in[0].parentid = 0;
      products->in[0].mech = 0;
      products->in[0].momentum = make_s_Momentum();
      products->in[0].momentum->px = evt.q1.Px();
      products->in[0].momentum->py = evt.q1.Py();
      products->in[0].momentum->pz = evt.q1.Pz();
      products->in[0].momentum->E = evt.q1.E();
      
      //PRODUCED PHOTON
      products->in[1].type = Gamma;
      products->in[1].pdgtype = 22;
      products->in[1].id = 2;
      products->in[1].parentid = 0;
      products->in[1].mech = 0;
      products->in[1].momentum = make_s_Momentum();
      products->in[1].momentum->px = evt.q2.Px();
      products->in[1].momentum->py = evt.q2.Py();
      products->in[1].momentum->pz = evt.q2.Pz();
      products->in[1].momentum->E = evt.q2.E(); 

      //PRODUCED PHOTON
      products->in[2].type = Gamma;
      products->in[2].pdgtype = 22;
      products->in[2].id = 3;
      products->in[2].parentid = 0;
      products->in[2].mech = 0;
      products->in[2].momentum = make_s_Momentum();
      products->in[2].momentum->px = evt.q3.Px();
      products->in[2].momentum->py = evt.q3.Py();
      products->in[2].momentum->pz = evt.q3.Pz();
      products->in[2].momentum->E = evt.q3.E();
      
      //PRODUCED PHOTON
      products->in[3].type = Gamma;
      products->in[3].pdgtype = 22;
      products->in[3].id = 4;
      products->in[3].parentid = 0;
      products->in[3].mech = 0;
      products->in[3].momentum = make_s_Momentum();
      products->in[3].momentum->px = evt.q4.Px();
      products->in[3].momentum->py = evt.q4.Py();
      products->in[3].momentum->pz = evt.q4.Pz();
      products->in[3].momentum->E = evt.q4.E();

      //PRODUCED PHOTON
      products->in[4].type = Gamma;
      products->in[4].pdgtype = 22;
      products->in[4].id = 5;
      products->in[4].parentid = 0;
      products->in[4].mech = 0;
      products->in[4].momentum = make_s_Momentum();
      products->in[4].momentum->px = evt.q5.Px();
      products->in[4].momentum->py = evt.q5.Py();
      products->in[4].momentum->pz = evt.q5.Pz();
      products->in[4].momentum->E = evt.q5.E();
      
      //PRODUCED PHOTON
      products->in[5].type = Gamma;
      products->in[5].pdgtype = 22;
      products->in[5].id = 6;
      products->in[5].parentid = 0;
      products->in[5].mech = 0;
      products->in[5].momentum = make_s_Momentum();
      products->in[5].momentum->px = evt.q6.Px();
      products->in[5].momentum->py = evt.q6.Py();
      products->in[5].momentum->pz = evt.q6.Pz();
      products->in[5].momentum->E = evt.q6.E();
     
      //PRODUCED Nucleus recoil
      //products->in[6].type = Helium;
      //products->in[6].pdgtype = 1000020040;
      if (evt.str_target == "He4" || evt.str_target == "Helium") {
	products->in[6].type = Helium;
	products->in[6].pdgtype = 1000020040;
      } else if (evt.str_target == "Be9" || evt.str_target == "Beryllium-9") {
	products->in[6].type = Be9;
	products->in[6].pdgtype = 1000040090;
      } else if (evt.str_target == "Proton") {
	products->in[6].type = Proton;
	products->in[6].pdgtype = 2212;
      } else if (evt.str_target == "Neutron") {
	products->in[6].type = Neutron;
	products->in[6].pdgtype = 2112;
      }
      products->in[6].id = 7;
      products->in[6].parentid = 0;
      products->in[6].mech = 0;
      products->in[6].momentum = make_s_Momentum();
      products->in[6].momentum->px = evt.q7.Px();
      products->in[6].momentum->py = evt.q7.Py();
      products->in[6].momentum->pz = evt.q7.Pz();
      products->in[6].momentum->E = evt.q7.E();
    }
    flush_s_HDDM(hddmEvt, ostream);

  }
};

#endif /* HDDMOUT_H_ */
