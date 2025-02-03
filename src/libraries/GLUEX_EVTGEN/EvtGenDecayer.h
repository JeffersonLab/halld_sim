
#ifndef __EVTGENDECAYER_H__
#define __EVTGENDECAYER_H__


#include <TLorentzVector.h>

#include "EvtGen/EvtGen.hh"

#include <vector>
#include <utility>

using namespace std;

// TODO: keep some global variable around that will check to see if this is initialized or not


class EvtGenDecayer {

	public:
		EvtGenDecayer(string in_user_decay_filename = "userDecay.dec");
		
		vector< pair< TLorentzVector, int > > decayParticle( const TLorentzVector& parent, int particleID );

	private:
		EvtGen *myGenerator = nullptr;
};



#endif   // __EVTGENDECAYER_H__
