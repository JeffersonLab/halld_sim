#if !defined(BREITWIGNERGENERATOR)
#define BREITWIGNERGENERATOR

#include <utility>
#include "TRandom.h"

using namespace std;

class BreitWignerGenerator
{
  
public:
  
  BreitWignerGenerator( int seed = 0 );
  BreitWignerGenerator( double mass, double width, int seed = 0 );
  
  void setSeed( unsigned int seed ){ m_randGen.SetSeed( seed ); }
  
  // output of the generation is a pair of doubles
  // the first is the mass and the second is the weight
  // to apply to this event to get back phase space
  //
  // values spanning the central fraction (optional argument)
  // of the distribution will be generated -- this allows
  // a mechanism to remove extreme values if desired
  pair< double, double > operator()( double fraction = 1 ) const;
  
  // returns the value of the PDF for some value of s
  double pdf( double s ) const;
  
private:
  
  mutable TRandom m_randGen;
    
  double random( double low, double hi ) const;
  
  static const double kPi;
  
  double m_mass;
  double m_width;
};

#endif

