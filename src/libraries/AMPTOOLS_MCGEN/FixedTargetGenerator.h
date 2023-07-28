#if !defined(FIXEDTARGETGENERATOR)
#define FIXEDTARGETGENERATOR

#include "TRandom.h"
#include "TLorentzVector.h"
#include "TH1.h"

#include "DecayChannelGenerator.h"
#include "BreitWignerGenerator.h"

#include <vector>

#define PI 3.14159

class Kinematics;

class FixedTargetGenerator {
  
public:
  
  FixedTargetGenerator(){}
  
  FixedTargetGenerator( double beamEnergy, double targetMass,
                        const vector< double >& uvMasses,
                        const vector< double >& lvMasses );
    
  Kinematics* generate() const;

  void setUpperVtxMasses( const vector< double >& uvMasses );
  void addUpperVertexBW( double mass, double width, double fraction );

  void setLowerVtxMasses( const vector< double >& lvMasses );
  void addLowerVertexBW( double mass, double width, double fraction );

  void setTargetMass( double targetMass );
  void setBeamEnergy( double eBeam );
  
  void setBeamP4( const TLorentzVector& p4 ){ m_beam = p4; }
  
  void addMomentumTransfer( double tSlope, double fraction );
  
private:

  void calculateLimits() const;

  pair< double, double > genUpperVertexMass() const;
  pair< double, double > genLowerVertexMass() const;
  
  vector< double > genPairMasses( double max, const vector< double >& masses ) const;

  vector< TLorentzVector > vertexGenP4( const TLorentzVector& vertexP4,
                                        const vector< double >& pairMasses,
                                        const vector< double >& masses ) const;
  
  double vertexPS( const vector< double >& pairMasses,
                   const vector< double >& masses ) const;

  pair< TLorentzVector, TLorentzVector >
       decay( const TLorentzVector& p4Initial, double m1, double m2 ) const;
  
  double pcm( double ecm, double m1, double m2 ) const;

  mutable TRandom m_randGen;
  
  TLorentzVector m_beam;
  TLorentzVector m_target;

  vector< double > m_uvMasses;
  vector< double > m_lvMasses;
  vector< double > m_tSlope;
  
  mutable double m_W;
  mutable double m_uvMax;
  mutable double m_uvMin;
  mutable double m_lvMax;
  mutable double m_lvMin;
  mutable double m_maxWeight;
  
  mutable bool m_limitsValid;
  
  DecayChannelGenerator m_upperDecayChannel;
  DecayChannelGenerator m_lowerDecayChannel;
  DecayChannelGenerator m_tSlopeChannel;
  
  vector< BreitWignerGenerator > m_upperBWGen;
  vector< BreitWignerGenerator > m_lowerBWGen;
  vector< double > m_tSlopeGen;
};

#endif
