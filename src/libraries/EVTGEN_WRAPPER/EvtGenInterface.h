//

#ifndef _EvtGenInterface_H_
#define _EvtGenInterface_H_

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


class EvtGenInterface 
{
	public:
		EvtGenInterface() {}
		~EvtGenInterface() {}
	
		void InitEvtGen();
		void Decay(hddm_s::HDDM * hddmevent);
		
		void SetPDGConversionMap(map<int,int> &in_map)  { PDGTYPE_CONVERSION_MAP = in_map; }
		void AddPDGConversionMap(int type1, int type2)  { PDGTYPE_CONVERSION_MAP[type1] = type2; }
		map<int,int> &GetPDGConversionMap(map<int,int> &in_map)  { return PDGTYPE_CONVERSION_MAP; }
		
		void SetUserDecay(string &infile) { USER_DECAY = infile; }
		string &GetUserDecay() { return USER_DECAY; }
		
		void SetGenSChannel(bool in_flag) { GEN_SCHANNEL = in_flag; }
		int GetGenSChannel(){ return GEN_SCHANNEL; }
		
	private:

		typedef struct {
			int id = -1;
			Particle_t type = Unknown;
			int pdgtype = -1;   
			bool decayed = false;
			TLorentzVector momentum;
			TVector3 vertex;
			hddm_s::ProductList::iterator hddmProduct;
		} gen_particle_info_t;


		void ParseVertices(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, int &max_particle_id);
		void DecayParticles(hddm_s::HDDM * hddmevent, vector< gen_particle_info_t > &particle_info, int &max_particle_id, int &vertex_id);

		EvtGen *myGenerator = nullptr;
		
		string USER_DECAY = "userDecay.dec";
		map<int,int> PDGTYPE_CONVERSION_MAP;
		
		// have to special-case some things for gen_schannel events
		bool GEN_SCHANNEL = false;
};

#endif // _EvtGenInterface_H_