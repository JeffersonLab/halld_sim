
#include "FixedTargetGenerator.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "Kinematics.h"

#include <iostream>

using namespace std;

FixedTargetGenerator::FixedTargetGenerator( double beamEnergy, double targetMass,
                                            const vector< double >& uvMasses,
                                            const vector< double >& lvMasses ) :
m_uvMasses( uvMasses ),
m_lvMasses( lvMasses ),
m_limitsValid( false )
{
  setBeamEnergy( beamEnergy );
  setTargetMass( targetMass );
  
  calculateLimits();
}

void
FixedTargetGenerator::addUpperVertexBW( double mass, double width, double fraction = 1 ){
  
  m_upperBWGen.push_back( BreitWignerGenerator( mass, width ) );
  m_upperDecayChannel.addChannel( m_upperBWGen.size()-1, fraction );
}

void
FixedTargetGenerator::addLowerVertexBW( double mass, double width, double fraction = 1 ){
  
  m_lowerBWGen.push_back( BreitWignerGenerator( mass, width ) );
  m_lowerDecayChannel.addChannel( m_lowerBWGen.size()-1, fraction );
}

void
FixedTargetGenerator::addMomentumTransfer( double tSlope, double fraction = 1 ){
  
  // we will always generate a negative argument:  exp( - t * tSlope  )
  // so take absolute value to avoid sign ambiguated
  m_tSlopeGen.push_back( fabs( tSlope ) );
  m_tSlopeChannel.addChannel( m_tSlopeGen.size(), fraction );
}

Kinematics*
FixedTargetGenerator::generate() const {

  if( !m_limitsValid ) calculateLimits();

  double lvM, uvM;
  vector< double > lvPairs;
  vector< double > uvPairs;

  double psWeight = 1;
  double seedWeight = 1;

  do {
    do {
      
      // first need to get two kinematically allowed
      // values for the upper and lower vertex masses
      do {
        
        seedWeight = 1;
        
        pair< double, double > tmp = genLowerVertexMass();
        lvM = tmp.first;
        seedWeight *= tmp.second;
        
        tmp = genUpperVertexMass();
        uvM = tmp.first;
        seedWeight *= tmp.second;
      }
      while( lvM + uvM > m_W );
      
      lvPairs = genPairMasses( m_lvMax, m_lvMasses );
      uvPairs = genPairMasses( m_uvMax, m_uvMasses );
    }
    
    // for single particle vertex, then the mass is guaranteed to be
    // within bounds -- otherwise compare the last pair + the last
    // stable particle with the total vertex mass
    while( ( ( m_lvMasses.size() > 1 ) &&
            ( lvPairs.back() + m_lvMasses.back() > lvM ) ) ||
           ( ( m_uvMasses.size() > 1 ) &&
            ( uvPairs.back() + m_uvMasses.back() > uvM ) ) );
        
    // put the total vertex mass on the end to make other
    // algorithms from here onward easier to write
    lvPairs.push_back( lvM );
    uvPairs.push_back( uvM );
    
    psWeight = pcm( m_W, lvM, uvM )/(4*PI*m_W);
    
    psWeight *= vertexPS( uvPairs, m_uvMasses );
    psWeight *= vertexPS( lvPairs, m_lvMasses );
    
    // there is a bug in the computation of the max weight
    // if this condition is not met
    assert( psWeight < m_maxWeight );
  }
  while( m_randGen.Uniform( 0, m_maxWeight ) > psWeight );
      
  // now generate the total upper and lower vertex four-vectors
  pair< TLorentzVector, TLorentzVector >
    top = decay( ( m_beam + m_target ), uvM, lvM );
  
  // from these generate the cascade of four-vectors from each vertex
  vector< TLorentzVector > uvP4 = vertexGenP4( top.first, uvPairs, m_uvMasses );
  vector< TLorentzVector > lvP4 = vertexGenP4( top.second, lvPairs, m_lvMasses );

  // concatenate the lists together
  vector< TLorentzVector > allP4 = uvP4;
  allP4.insert( allP4.end(), lvP4.begin(), lvP4.end() );
    
  return new Kinematics( allP4, seedWeight );
}

void
FixedTargetGenerator::setBeamEnergy( double energy ){
  
  m_beam = TLorentzVector( 0, 0, energy, energy );
  m_limitsValid = false;
}

void
FixedTargetGenerator::setTargetMass( double mass ){
  
  m_target = TLorentzVector( 0, 0, 0, mass );
  m_limitsValid = false;
}


// the second value in the pair is the weight needed
// to restore a uniform distribution in mass
pair< double, double >
FixedTargetGenerator::genUpperVertexMass() const {
  
  // single particle vertex -- return the mass with
  // a weight of 1
  if( m_uvMasses.size() == 1 )
    return pair< double, double >( m_uvMasses[0], 1 );
  
  // no BW's to seed the generation:  return uniform
  if( m_upperBWGen.size() == 0 )
    return pair< double, double >( m_randGen.Uniform( m_uvMin, m_uvMax ), 1 );

  pair< double, double > massWeight;
  do{

    massWeight =
      m_upperBWGen[m_upperDecayChannel()]();
    
    // the BW generator returns values that are flat
    // in 2-body phase space (s) reweight to make it
    // them flat in mass by multiplying by 1/M
    
    massWeight.second /= massWeight.first;
  }
  while( ( massWeight.first < m_uvMin ) ||
         ( massWeight.first > m_uvMax ) );
 
  return massWeight;
}

// the second value in the pair is the weight needed
// to restore a uniform distribution in mass
pair< double, double >
FixedTargetGenerator::genLowerVertexMass() const {
  
  // single particle vertex -- return the mass with
  // a weight of 1
  if( m_lvMasses.size() == 1 )
    return pair< double, double >( m_lvMasses[0], 1 );
  
  // no BW's to seed the generation:  return uniform
  if( m_lowerBWGen.size() == 0 )
    return pair< double, double >( m_randGen.Uniform( m_lvMin, m_lvMax ), 1 );
  
  pair< double, double > massWeight;
  do{

    massWeight =
      m_lowerBWGen[m_lowerDecayChannel()]();
    
    // the BW generator returns values that are flat
    // in 2-body phase space (s) reweight to make it
    // them flat in mass by multiplying by 1/M
    
    massWeight.second /= massWeight.first;
  }
  while( ( massWeight.first < m_lvMin ) ||
         ( massWeight.first > m_lvMax ) );
 
  return massWeight;
}

void
FixedTargetGenerator::calculateLimits() const {
  
  // this is the sqrt(s):
  m_W = ( m_beam + m_target ).M();
  
  double lvSum = 0;
  for( unsigned int i = 0; i < m_lvMasses.size(); ++i ) lvSum += m_lvMasses[i];

  double uvSum = 0;
  for( unsigned int i = 0; i < m_uvMasses.size(); ++i ) uvSum += m_uvMasses[i];
  
  m_lvMin = lvSum;
  m_lvMax = m_W - uvSum;
  
  m_uvMin = uvSum;
  m_uvMax = m_W - lvSum;
  
  // now calculate the max weight by computing the products of the maximum
  // relative momentum at every step of the decay
  
  m_maxWeight = pcm( m_W, m_lvMin, m_uvMin )/(4*PI*m_W);
  
  // for each of the vertices, we are going to step backwards starting
  // from the final two stable particles -- the maximum momentum for these
  // two occurs when all other particles "above" them are produced at
  // rest
  
  double eMax = m_uvMax - uvSum + m_uvMasses[0];
  double minMass = 0;
  for( unsigned int i = 1; i < m_uvMasses.size(); ++i ){
  
    minMass += m_uvMasses[i-1];
    eMax += m_uvMasses[i];
    
    m_maxWeight *= pcm( eMax, minMass, m_uvMasses[i] );
  }

  eMax = m_lvMax - lvSum + m_lvMasses[0];
  minMass = 0;
  for( unsigned int i = 1; i < m_lvMasses.size(); ++i ){
  
    minMass += m_lvMasses[i-1];
    eMax += m_lvMasses[i];
    
    m_maxWeight *= pcm( eMax, minMass, m_lvMasses[i] );
  }
  
  m_limitsValid = true;
}

vector< double >
FixedTargetGenerator::genPairMasses( double max, const vector< double >& masses ) const {
  
  if( masses.size() == 1 ) return masses;
  
  // this creates a vector called pairMasses that holds the invariant
  // masses of the intermediate states -- one should have the decays:
  //    pairMasses[i] -> pairMasses[i-1] + masses[i]
  // for convenience we set pairMasses[0] = masses[0]
  //
  
  double min = 0;
  for( unsigned int i = 0; i < masses.size(); ++i ) min += masses[i];
  
  // the invaraint masses at the top level have already been generated
  // this leaves N - 2 masses to pick, it is helpful to load the first
  // element of pairMasses with the first stable paricle mass
  vector< double > pairMasses( masses.size() - 1 );

  // the method used below for efficiently generating the mass pairs
  // is called the Raubold-Lynch method and it is extensively
  // documented in CERN technical note by Fred James

  vector< double > r( masses.size() - 1 );
  for( unsigned int i = 1; i < r.size(); ++i ) r[i] = m_randGen.Uniform();
  // the first element sorted r will be zero
  sort( r.begin(), r.end() );

  pairMasses[0] = masses[0];

  double thisMin = masses[0];
  double thisDelta = max - min;
  
  for( unsigned int i = 1; i < r.size(); ++i ){

    thisMin += masses[i];
    thisDelta += masses[i];

    pairMasses[i] = r[i]*thisDelta + thisMin;
  }

  return pairMasses;
}

vector< TLorentzVector >
FixedTargetGenerator::vertexGenP4( const TLorentzVector& vertexP4,
                                   const vector< double >& pairMasses,
                                   const vector< double >& masses ) const {
  
  if( masses.size() == 1 )
    return vector< TLorentzVector >( 1, vertexP4 );

  pair< TLorentzVector, TLorentzVector > aPair;
  vector< TLorentzVector > p4List( masses.size() );
  
  // traverse from the back of the pairMass list, the last element
  // should have the total mass of the vertex:  vertexP4.M() so
  // for generating decays we need to start with the 2nd to last
  // element -- at the other end note pairMasses[0] == masses[0]

  TLorentzVector parentP4 = vertexP4;
  
  for( int i = pairMasses.size() - 2; i >= 0; --i ){
    
    aPair = decay( parentP4, pairMasses.at(i), masses.at(i+1) );

    // we save the second particle keeping the same ordering
    // as in the array of masses and decay the first
    p4List[i+1] = aPair.second;
    parentP4 = aPair.first;
  }
  
  // now save the final stable particle
  p4List[0] = aPair.first;

  return p4List;
}

double
FixedTargetGenerator::vertexPS( const vector< double >& pairMasses,
                                const vector< double >& masses ) const {

  if( masses.size() == 1 ) return 1;
    
  double psSize = 1;
  for( unsigned int i = 1; i < pairMasses.size(); ++i ){
    
    // note that pairMasses[0] == masses[0]
    psSize *= pcm( pairMasses.at(i), pairMasses.at(i-1), masses.at(i) );
  }

  return psSize;
}

double
FixedTargetGenerator::pcm( double ecm, double m1, double m2 ) const {
 
  // the Kallen function
  double lambda = ( ecm*ecm - m1*m1 - m2*m2 )*
                  ( ecm*ecm - m1*m1 - m2*m2 ) -
                  4*m1*m1*m2*m2;
  return sqrt( lambda ) / ( 2 * ecm );
}

pair< TLorentzVector, TLorentzVector >
FixedTargetGenerator::decay( const TLorentzVector& p4Initial,
                             double m1, double m2 ) const {
  
  pair< TLorentzVector, TLorentzVector > p4Final;
  
  // first generate the decay of back to back particles
  // in the center of momentum frame
  p4Final.first.SetXYZM( 0, 0, pcm( p4Initial.M(), m1, m2 ), m1 );
  p4Final.first.RotateY( acos( m_randGen.Uniform( -1, 1 ) ) );
  p4Final.first.RotateZ( m_randGen.Uniform( 0, 2*PI ) );
  p4Final.second.SetVectM( -p4Final.first.Vect(), m2 );
  
  // then boost to the frame of initial state
  p4Final.first.Boost( p4Initial.BoostVector() );
  p4Final.second.Boost( p4Initial.BoostVector() );
  
  return p4Final;
}

