#if !defined(FIXEDTARGETGENERATOR)
#define FIXEDTARGETGENERATOR

#include "TRandom.h"
#include "TLorentzVector.h"
#include "TH1.h"

#include "AMPTOOLS_MCGEN/DecayChannelGenerator.h"
#include "AMPTOOLS_MCGEN/BreitWignerGenerator.h"

#include <vector>

#define PI 3.14159

class Kinematics;

class FixedTargetGenerator {
  
public:
  
  enum { kMomentumTransfer = 1, kUpperVtxMass = 2, kLowerVtxMass = 4 };
  
  FixedTargetGenerator();
  
  // this constructor assumes the beam is a photon (go GlueX!)  --
  // more general cases can be handled by setBeamP4 below
  FixedTargetGenerator( double photonBeamEnergy, double targetMass,
                        const vector< double >& uvMasses,
                        const vector< double >& lvMasses );

  void setSeed( unsigned int seed ){ m_randGen.SetSeed( seed ); }
  
  Kinematics* generate() const;

  // These functions setup the inital state.
  void setBeamEnergy( double eBeam );
  void setBeamP4( const TLorentzVector& p4 );
  void setTargetMass( double targetMass );

  // Use to set the masses of the stable particles produced
  // at the upper and lower verticies
  void setUpperVtxMasses( const vector< double >& uvMasses );
  void setLowerVtxMasses( const vector< double >& lvMasses );

  // This is used to seed generation with a BW to allow the user
  // to concentrante events in regions of mass where the amplitude
  // might peak.  If the kUpperVtxMass and/or kLowerVtxMass bits
  // are turned on in the reweight mask then weights will be provided
  // to smooth out these peaks in one or both distributions.
  void addUpperVtxBW( double mass, double width, double fraction );
  void addLowerVtxBW( double mass, double width, double fraction );

  // specify a range of mass to generate -- otherwise the full
  // kinematic limits are used
  void setUpperVtxRange( double min, double max );
  void setLowerVtxRange( double min, double max );
  
  // these should be absolute values:  |t|
  void setMomentumTransferRange( double min, double max );
  void addMomentumTransfer( double tSlope, double fraction );

  // The mask a bitmask formed with the enums above -- the
  // default behavior is achieved by:
  //
  // setReweightMask( kMomentumTransfer | kUpperVtxMass | kLowerVtxMass );
  //
  // In this case if the user inputs either BW distributions for
  // the vertices or exponential t slopes then then the events will
  // be generated with a weight that recovers pure phase space.
  //
  // A desirable option for some applications may be the following:
  //
  // addMomentumTransfer( 2.0 );
  // setReweightMask( kUpperVtxMass | kLowerVtxMass );
  //
  // In this case events will be generated with e^-2*|t| and, additional
  // upper vertex and lower vertex BW's as desired.  The weight with the
  // events will unweight the BW's but not the t-distribution.  The
  // corresponding mass distributions then will no longer be pure phase
  // space but what one gets when phase space is weighted by e^-2*|t|
  
  void setReweightMask( unsigned int mask );

private:

  void calculateLimits() const;

  double pcm( double ecm, double m1, double m2 ) const;

  pair< double, double > genUpperVertexMass() const;
  pair< double, double > genLowerVertexMass() const;
  pair< double, double > genCosThetaCM( double uvM, double lvM, double* t ) const;
  
  vector< double > genPairMasses( double max, const vector< double >& masses ) const;

  vector< TLorentzVector > vertexGenP4( const TLorentzVector& vertexP4,
                                        const vector< double >& pairMasses,
                                        const vector< double >& masses ) const;
  
  double vertexPS( const vector< double >& pairMasses,
                   const vector< double >& masses ) const;

  pair< TLorentzVector, TLorentzVector >
       decay( const TLorentzVector& p4Initial, double m1, double m2 ) const;
  pair< TLorentzVector, TLorentzVector >
       decay( const TLorentzVector& p4Initial, double m1, double m2, double cosThetaCM ) const;
  
  mutable TRandom m_randGen;
  
  TLorentzVector m_beam;
  TLorentzVector m_target;

  vector< double > m_uvMasses;
  vector< double > m_lvMasses;
  vector< double > m_tSlope;
  
  double m_lvMinUser;
  double m_lvMaxUser;
  double m_uvMinUser;
  double m_uvMaxUser;
  double m_tMagMinUser;
  double m_tMagMaxUser;
  unsigned int m_reweightMask;
  
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
  vector< double > m_tSlopeVec;
};

#endif
