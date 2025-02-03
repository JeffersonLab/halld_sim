
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"

#include "AMPTOOLS_MCGEN/FixedTargetGenerator.h"
#include "IUAmpTools/Kinematics.h"

#include <iostream>

using namespace std;

FixedTargetGenerator::FixedTargetGenerator( int seed ) :
m_seed( seed ),
m_lvMinUser( 0 ),
m_lvMaxUser( 1E9 ),
m_uvMinUser( 0 ),
m_uvMaxUser( 1E9 ),
m_tMagMinUser( 0 ),
m_tMagMaxUser( 1E9 ),
// by default we reweight any seed distributions put in by the user
// to recover phase space
m_reweightMask( kUpperVtxMass | kLowerVtxMass | kMomentumTransfer ),
m_limitsValid( false )
{
  setSeed( seed );
}

FixedTargetGenerator::FixedTargetGenerator( double photonBeamEnergy, double targetMass,
                                            const vector< double >& uvMasses,
                                            const vector< double >& lvMasses,
                                            int seed ) :
m_seed( seed ),
m_uvMasses( uvMasses ),
m_lvMasses( lvMasses ),
m_lvMinUser( 0 ),
m_lvMaxUser( 1E9 ),
m_uvMinUser( 0 ),
m_uvMaxUser( 1E9 ),
m_tMagMinUser( 0 ),
m_tMagMaxUser( 1E9 ),
// by default we reweight any seed distributions put in by the user
// to recover phase space
m_reweightMask( kUpperVtxMass | kLowerVtxMass | kMomentumTransfer ),
m_limitsValid( false )
{

  setBeamEnergy( photonBeamEnergy );
  setTargetMass( targetMass );
  setSeed( seed );
  
  calculateLimits();
}

void
FixedTargetGenerator::setSeed( int seed ){
  
  // save the seed in case we create more generators
  // that need to be seeded
  m_seed = seed;
  
  // now seed every generator that this class is
  // is responsible for
  m_randGen.SetSeed( seed );

  m_upperDecayChannel.setSeed( seed );
  m_lowerDecayChannel.setSeed( seed );
  m_tSlopeChannel.setSeed( seed );
  
  for( vector< BreitWignerGenerator >::iterator bwGen =
      m_upperBWGen.begin(); bwGen != m_upperBWGen.end();
      ++bwGen ){
    
    (*bwGen).setSeed( seed );
  }

  for( vector< BreitWignerGenerator >::iterator bwGen =
      m_lowerBWGen.begin(); bwGen != m_lowerBWGen.end();
      ++bwGen ){
    
    (*bwGen).setSeed( seed );
  }
}

void
FixedTargetGenerator::setUpperVtxMasses( const vector< double >& uvMasses ){
  
  m_uvMasses = uvMasses;
  m_limitsValid = false;
}

void
FixedTargetGenerator::setLowerVtxMasses( const vector< double >& lvMasses ){
  
  m_lvMasses = lvMasses;
  m_limitsValid = false;
}

void
FixedTargetGenerator::addUpperVtxBW( double mass, double width, double fraction ){
  
  m_upperBWGen.push_back( BreitWignerGenerator( mass, width, m_seed ) );
  m_upperDecayChannel.addChannel( m_upperBWGen.size()-1, fraction );
}

void
FixedTargetGenerator::addLowerVtxBW( double mass, double width, double fraction ){
  
  m_lowerBWGen.push_back( BreitWignerGenerator( mass, width, m_seed ) );
  m_lowerDecayChannel.addChannel( m_lowerBWGen.size()-1, fraction );
}

void
FixedTargetGenerator::addMomentumTransfer( double tSlope, double fraction ){
  
  // we will always generate a negative argument:  exp( - t * tSlope  )
  // so take absolute value to avoid sign ambiguated

  m_tSlopeVec.push_back( fabs( tSlope ) );
  m_tSlopeChannel.addChannel( m_tSlopeVec.size()-1, fraction );
}

void
FixedTargetGenerator::setBeamEnergy( double energy ){
  
  m_beam = TLorentzVector( 0, 0, energy, energy );
  m_limitsValid = false;
}

void
FixedTargetGenerator::setBeamP4( const TLorentzVector& p4 ){
  
  // if the beam is not aligned with z, then the logic in the
  // algorithm for generating cos(theta) for user-specified
  // momentum transfer will not work... it isn't general enough
  // to handle those cases
  if( p4.Perp() > 0 ){
    
    cerr << "FixedTargetGenerator: only beams in the z direction are supported" << endl;
    exit( 1 );
  }
  
  m_beam = p4;
}

void
FixedTargetGenerator::setTargetMass( double mass ){
  
  m_target = TLorentzVector( 0, 0, 0, mass );
  m_limitsValid = false;
}

void
FixedTargetGenerator::setUpperVtxRange( double min, double max ){
  
  m_uvMinUser = min;
  m_uvMaxUser = max;
  m_limitsValid = false;
}

void
FixedTargetGenerator::setLowerVtxRange( double min, double max ){
  
  m_lvMinUser = min;
  m_lvMaxUser = max;
  m_limitsValid = false;
}

void
FixedTargetGenerator::setMomentumTransferRange( double min, double max ){
  
  // try to catch confusion with t and |t|...
  //   a more informative message might be nice
  assert( fabs( max ) > fabs( min ) );
  
  m_tMagMinUser = fabs( min );
  m_tMagMaxUser = fabs( max );
}

void
FixedTargetGenerator::setReweightMask( unsigned int mask ){
  
  // we are only using the 3 lowest bits of the int
  assert( mask < 8 );
  m_reweightMask = mask;
}

// **
// const member functions:

Kinematics*
FixedTargetGenerator::generate( bool includeBeam ) const {

  if( !m_limitsValid ) calculateLimits();

  double lvM, uvM, cosThetaCM;
  vector< double > lvPairs;
  vector< double > uvPairs;

  double psWeight = 1;
  double seedWeight = 1;

  do {
    do {
      
      // first need to get two kinematically allowed
      // values for the upper and lower vertex masses
      // and the scattering angle in the CM frame
      
      double t = 0;
      do {
        
        seedWeight = 1;
        
        pair< double, double > tmp = genLowerVertexMass();
        lvM = tmp.first;
        if( m_reweightMask & kLowerVtxMass ) seedWeight *= tmp.second;
        
        tmp = genUpperVertexMass();
        uvM = tmp.first;
        if( m_reweightMask & kUpperVtxMass ) seedWeight *= tmp.second;
      
        tmp = genCosThetaCM( uvM, lvM, &t );
        cosThetaCM = tmp.first;
        if( m_reweightMask & kMomentumTransfer ) seedWeight *= tmp.second;
      }
      while( ( lvM + uvM > m_W ) ||
             ( fabs( cosThetaCM ) > 1  ) );
      
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
    
    // here we are computing the weights to recover N-body
    // phase space from events that are drawn uniform in
    // mass at the upper and lower vertex -- this means the
    // while( ) statement at the bottom will do accept/reject
    // on the phase space weight only.  This produces events
    // of unity weight distributed by phase space.  These events
    // may be further modified by user "seeds" and weighted
    // by seedWeight to recover phase space.  True phase space
    // will only be recovered if m_reweightMask = 7 (all bits on).
    
    psWeight = pcm( m_W, lvM, uvM )/(4*TMath::Pi()*m_W);
    
    psWeight *= vertexPS( uvPairs, m_uvMasses );
    psWeight *= vertexPS( lvPairs, m_lvMasses );
    
    // there is a bug in the computation of the max weight
    // if this condition is not met
    assert( psWeight <= m_maxWeight );
  }
  while( m_randGen.Uniform( 0, m_maxWeight ) > psWeight );
      
  // now generate the total upper and lower vertex four-vectors
  pair< TLorentzVector, TLorentzVector >
    top = decay( ( m_beam + m_target ), uvM, lvM, cosThetaCM );
    
  // from these generate the cascade of four-vectors from each vertex
  vector< TLorentzVector > uvP4 = vertexGenP4( top.first, uvPairs, m_uvMasses );
  vector< TLorentzVector > lvP4 = vertexGenP4( top.second, lvPairs, m_lvMasses );

  // concatenate the lists together
  vector< TLorentzVector > allP4;
  if( includeBeam ) allP4.push_back( m_beam );
  allP4.insert( allP4.end(), lvP4.begin(), lvP4.end() );
  allP4.insert( allP4.end(), uvP4.begin(), uvP4.end() );
    
  return new Kinematics( allP4, seedWeight );
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
    
    massWeight = m_upperBWGen[m_upperDecayChannel()]();
    
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

    massWeight = m_lowerBWGen[m_lowerDecayChannel()]();
    
    // the BW generator returns values that are flat
    // in 2-body phase space (s) reweight to make it
    // them flat in mass by multiplying by 1/M
    
    massWeight.second /= massWeight.first;
  }
  while( ( massWeight.first < m_lvMin ) ||
         ( massWeight.first > m_lvMax ) );
 
  return massWeight;
}

pair< double, double >
FixedTargetGenerator::genCosThetaCM( double uvMass, double lvMass, double* t ) const {

  pair< double, double > cosThetaWeight;

  double eBeamCM = ( m_W*m_W + m_beam.M2() - m_target.M2() ) / ( 2*m_W );
  double puvCM = pcm( m_W, uvMass, lvMass );
  double euvCM = sqrt( puvCM*puvCM + uvMass*uvMass );
  
  do{
    
    if( m_tSlopeVec.size() == 0 ) {
      
      // this is pure phase space:
      cosThetaWeight.first = m_randGen.Uniform( -1, 1 );
      cosThetaWeight.second = 1;
      
      // no seeding with exponentials -- compute t from cos( theta )
      *t = cosThetaWeight.first * ( 2 * pcm( m_W, m_beam.M(), m_target.M() ) * puvCM ) -
      ( 2*eBeamCM*euvCM - m_beam.M2() - uvMass*uvMass );
    }
    else{
      
      // if we get here, we're going to generate t and then use it compute
      // cos( theta ) prior to returning
      
      // get a tSlope according to the mix of slopes specified by the user
      double tSlope = m_tSlopeVec[m_tSlopeChannel()];
      
      // the ROOT generator generates exponential decay so tSlope should be positive
      *t = -m_randGen.Exp( 1./tSlope );
      
      // ref: Martin and Spearman Eq. 4.81b
      cosThetaWeight.first = ( *t + 2*eBeamCM*euvCM - m_beam.M2() - uvMass*uvMass ) /
      ( 2 * pcm( m_W, m_beam.M(), m_target.M() ) * puvCM );
      
      // dN / d(cos theta ) = 2 q q' [ dN / dt ] ... to make the LHS
      // uniform (phase space) we weight by the inverse of the RHS
      // (t is negative in this expression and tSlope is positive)
      cosThetaWeight.second = 1 /
         ( 2 * pcm( m_W, m_beam.M(), m_target.M() ) * puvCM * exp( *t * tSlope ) );
    }
  }
  while( ( fabs( *t ) < m_tMagMinUser ) ||
         ( fabs( *t ) > m_tMagMaxUser ) );
    
  return cosThetaWeight;
}


void
FixedTargetGenerator::calculateLimits() const {
  
  // this is the sqrt(s):
  m_W = ( m_beam + m_target ).M();
  
  double lvSum = 0;
  for( unsigned int i = 0; i < m_lvMasses.size(); ++i ) lvSum += m_lvMasses[i];

  double uvSum = 0;
  for( unsigned int i = 0; i < m_uvMasses.size(); ++i ) uvSum += m_uvMasses[i];
  
  assert( m_W >= lvSum + uvSum && "Sum of particle masses can not be larger than center of mass energy");
  // set the generation limits considering also any user-specified
  // range of masses that may be more restrictive than the kinematic limits
  
  m_lvMin = ( m_lvMinUser > lvSum ? m_lvMinUser : lvSum );
  m_lvMax = ( m_lvMaxUser < m_W - uvSum ? m_lvMaxUser : m_W - uvSum );
  
  m_uvMin = ( m_uvMinUser > uvSum ? m_uvMinUser : uvSum );
  m_uvMax = ( m_uvMaxUser < m_W - lvSum ? m_uvMaxUser : m_W - lvSum );
  
  // now calculate the max weight by computing the products of the maximum
  // relative momentum at every step of the decay
  
  m_maxWeight = pcm( m_W, m_lvMin, m_uvMin )/(4*TMath::Pi()*m_W);

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
FixedTargetGenerator::genPairMasses( double max,
                                     const vector< double >& masses ) const {
  
  if( masses.size() == 1 ) return masses;
  
  // this creates a vector called pairMasses that holds the invariant
  // masses of the intermediate states -- one should have the decays:
  //    pairMasses[i] -> pairMasses[i-1] + masses[i]
  // for convenience we set pairMasses[0] = masses[0]
  
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
  
  return decay( p4Initial, m1, m2, m_randGen.Uniform( -1, 1 ) );
}

pair< TLorentzVector, TLorentzVector >
FixedTargetGenerator::decay( const TLorentzVector& p4Initial,
                             double m1, double m2, double cosTheta ) const {
  
  // this method allows specification of cosTheta in the convention
  // where theta is angle w.r.t. z-axis in the rest frame of
  // p4Initial... doing so is probably only generally useful when
  // p4Initial is along z (in this case theta is the more-standard
  // helicity angle) -- use caution if you decide to copy this
  // function for genreal use elsewhere!!
  
  pair< TLorentzVector, TLorentzVector > p4Final;
  
  // first generate the decay of back to back particles
  // in the center of momentum frame
  p4Final.first.SetXYZM( 0, 0, pcm( p4Initial.M(), m1, m2 ), m1 );

  // then rotate appropriately
  p4Final.first.RotateY( acos( cosTheta ) );
  p4Final.first.RotateZ( m_randGen.Uniform( 0, 2*TMath::Pi() ) );
  p4Final.second.SetVectM( -p4Final.first.Vect(), m2 );
  
  // then boost to the frame of initial state
  p4Final.first.Boost( p4Initial.BoostVector() );
  p4Final.second.Boost( p4Initial.BoostVector() );
  
  return p4Final;
}

