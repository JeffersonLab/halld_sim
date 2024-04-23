#if !defined(DECAYCHANNELGENERATOR)
#define DECAYCHANNELGENERATOR

#include <vector>
#include "TRandom.h"

using namespace std;

class DecayChannelGenerator {
  
public:
  
  DecayChannelGenerator( int seed = 0 );
  
  void addChannel( unsigned int channelNum, double bf );
  
  void setSeed( unsigned int seed ){ m_randGen.SetSeed( seed ); }

  unsigned int operator()() const;
  const vector< unsigned int >& availableChannels() const { return m_index; }
  double getProb( unsigned int channelNum ) const;
  
private:
  
  mutable TRandom m_randGen;
  
  mutable double m_bfTotal;
  mutable bool m_probRenormalized;
  
  mutable vector< double > m_upperBound;
  mutable vector< double > m_prob;
  
  vector< unsigned int > m_index;
};

#endif
