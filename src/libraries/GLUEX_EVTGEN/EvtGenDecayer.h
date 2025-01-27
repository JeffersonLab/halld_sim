
#ifndef __EVTGENDECAYER_H__
#define __EVTGENDECAYER_H__


#include <TLorentzVector.h>

#include "EvtGen/EvtGen.hh"

#include <vector>
#include <utility>

using namespace std;


class EvtGenDecayer {

	public:
		EvtGenDecayer();
		
		vector< pair< TLorentzVector, int > > decayParticle( const TLorentzVector& parent, int particleID );

	private:
		EvtGen *myGenerator = nullptr;
};



#endif   // __EVTGENDECAYER_H__
