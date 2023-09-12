#if !defined(DECAYCHANNELGENERATOR)
#define DECAYCHANNELGENERATOR

#include <vector>

using namespace std;

class DecayChannelGenerator {
  
public:
  
  DecayChannelGenerator();
  
  void addChannel( unsigned int channelNum, double bf );
  
  unsigned int operator()() const;
  const vector< unsigned int >& availableChannels() const { return m_index; }
  double getProb( unsigned int channelNum ) const;
  
private:
  
  mutable double m_bfTotal;
  mutable bool m_probRenormalized;
  
  mutable vector< double > m_upperBound;
  mutable vector< double > m_prob;
  vector< unsigned int > m_index;
};

#endif
