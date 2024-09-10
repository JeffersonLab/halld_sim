
#include <iostream>
#include <stdlib.h>

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "particleType.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

const double ProductionMechanism::kPi = 3.14159;

using namespace std;

ProductionMechanism::ProductionMechanism( Recoil recoil, Type type, double slope, int seed ) :
m_type( type ),
m_lowMass( 0 ),
m_highMass( 0 ),
m_slope( slope ),
m_lowT( 0 ),
m_highT( 12 ),
m_lastWeight( 1. )
{	
  kMproton=ParticleMass(Proton);
  kMneutron=ParticleMass(Neutron);
  kMlambda=ParticleMass(Lambda);
  kMsigmaPlus=ParticleMass(SigmaPlus);
  // kMZ = 108.;      //  mass of Sn116 
  kMZ = 208.*0.931494;      //  use mass of Pb as it is in the particle table
  kMPion = ParticleMass(PiPlus);
  kMPi0 = ParticleMass(Pi0);
  kMKaon = ParticleMass(KPlus);
  
  isBaryonResonance = false;

  // initialize pseudo-random generator
  gRandom->SetSeed(seed);

  switch( recoil ){
    // I'm sure the distinction between these doesn't matter!  
  case kProton:  m_recMass = kMproton; break; //old value: 0.9382
  case kNeutron: m_recMass = kMneutron; break; //old value: 0.9395
  case kLambda: m_recMass = kMlambda; break;
  case kSigmaPlus: m_recMass = kMsigmaPlus; break;
  case kZ: m_recMass = kMZ; break; //default to Sn116/Pb
  case kPion: m_recMass = kMPion; isBaryonResonance = true; break;
  case kPi0: m_recMass = kMPi0; isBaryonResonance = true; break;
  case kKaon: m_recMass = kMKaon; isBaryonResonance = true; break;
  default:       m_recMass = kMproton; break; //old value: 0.9382
  }
}

void
ProductionMechanism::setMassRange( double low, double high ){
	
	m_lowMass = low;
	m_highMass = high;
}

void
ProductionMechanism::setTRange( double low, double high ){
	
	m_lowT = low;
	m_highT = high;
}

void 
ProductionMechanism::setGeneratorType( Type type ){
  
  m_type = type;
}

void
ProductionMechanism::setRecoilMass( double recMass ){

	m_recMass = recMass;
}

TLorentzVector
ProductionMechanism::produceResonance( const TLorentzVector& beam ){

  TLorentzVector target( 0, 0, 0, kMproton );
  
  TLorentzRotation lab2cmBoost( -( target + beam ).BoostVector() );
  TLorentzRotation cm2labBoost( ( target + beam ).BoostVector() );
  
  double cmEnergy = ( lab2cmBoost * ( target + beam ) ).E();
  double beamMomCM = cmMomentum( cmEnergy, beam.M(), target.M() );
  
  // double exptMax = exp(-1.)/m_slope;   // Elton 8/19/2016.  t*exp(Bt)
  double exptMax = 1;   // remove factor of t for rho production (no spin flip). set this value for exp(Bt)
  
  double t, tMin, tMax, resMass, resMomCM;

  // First generate flat mass distribution
  do // the resonance mass cannot be larger than CM energy - recoil mass
    resMass = generateMass();
  while ( cmEnergy < resMass + m_recMass );

  // Then generate t accordingly
  do {
    resMomCM  = cmMomentum( cmEnergy, resMass, m_recMass );
    
    tMin = 0;
    tMax = 4. * beamMomCM * resMomCM;
    
    double tlow(tMin), thigh(tMax);
    if ( tMin < m_lowT ) tlow=m_lowT;
    if ( m_highT < tMax) thigh=m_highT;
    t = random( tlow, thigh ); 
  } 
  // while( random( 0., exptMax ) > t*exp(-m_slope*t) );   // Elton 8/19/2016.  t*exp(Bt)
  while( random( 0., exptMax ) > exp(-m_slope*t) );   // remove factor of t for rho production (no spin flip). Set this line for exp(Bt)
  
  TVector3 resonanceMomCM;
  if(isBaryonResonance){
	resonanceMomCM.SetMagThetaPhi( resMomCM,
					kPi-acos( 1. - 2.*t/tMax ), // opposite of what it would be for meson resonances
					random( -kPi, kPi ) );
  }
  else{
	resonanceMomCM.SetMagThetaPhi( resMomCM,
					acos( 1. - 2.*t/tMax ),
					random( -kPi, kPi ) );
  }
  
  TLorentzVector resonanceCM( resonanceMomCM, 
			      sqrt( resonanceMomCM.Mag2() +
                                    resMass * resMass ) );
  
  return cm2labBoost * resonanceCM;
}
TLorentzVector
ProductionMechanism::produceResonanceZ ( const TLorentzVector& beam){
  /* This method is modeled after produceResonance, which assumes a proton target and exponential t dependence
     This method is intended for use with a high Z target in Primakoff production.  Elton 4/14/2017

   */
	
	TLorentzVector target( 0, 0, 0, kMZ);
	
	TLorentzRotation lab2cmBoost( -( target + beam ).BoostVector() );
	TLorentzRotation cm2labBoost( ( target + beam ).BoostVector() );
	
	double cmEnergy = ( lab2cmBoost * ( target + beam ) ).E();
	double beamMomCM = cmMomentum( cmEnergy, beam.M(), target.M() );

	// double exptMax = exp(-1.)/m_slope;   // Elton 8/19/2016.  t*exp(Bt)
        double exptMax = 1;   // remove factor of t for rho production (no spin flip). set this value for exp(Bt)

	double t, tMaxkin, tMax, resMass, resMomCM;
	// generate the t-distribution. t is positive here (i.e. should be -t)


  do {
    resMass = generateMass();
    resMomCM  = cmMomentum( cmEnergy, resMass, m_recMass );
  
    tMaxkin = 4. * beamMomCM * resMomCM;
    tMax = 0.2;   // restrict max to make more efficient for Primakoff generation (about 2. deg at 0.05 GeV-2)
    // tMax = 1.;   // restrict max to make more efficient for Primakoff generation
    t = random( 0, tMax ); 
  } 
  // while( random( 0., exptMax ) > t*exp(-m_slope*t) );   // Elton 8/19/2016.  t*exp(Bt)
  while( random( 0., exptMax ) > exp(-m_slope*t) );   // remove factor of t for rho production (no spin flip). Set this line for exp(Bt)

  // cout << endl << "produceResonanceZ, resMomCM=" << resMomCM << " resMass=" << resMass << " t=" << t << " tMax=" << tMax << " cmEnergy=" << cmEnergy << " kMZ=" << kMZ << endl;

	TVector3 resonanceMomCM;
	double thetaCM = 2.*sqrt(t/tMaxkin); // acos( 1. - 2.*t/tMax ) -> use small angle approximation to avoid roundoff.
	// double thetaCM = acos( 1. - 2.*t/tMaxkin );
	double phiCM = random( -kPi, kPi ); 

	resonanceMomCM.SetMagThetaPhi( resMomCM, thetaCM, phiCM);
	
	TLorentzVector resonanceCM( resonanceMomCM, 
                               sqrt( resonanceMomCM.Mag2() +
                                    resMass * resMass ) );
	// resonanceCM.Print();
	
	return cm2labBoost * resonanceCM;
}

void 
ProductionMechanism::addResonance( double mass, double width, double crossSec ){
  
  m_decGen.addChannel( m_bwGen.size(), crossSec );
  m_bwGen.push_back( BreitWignerGenerator( mass, width ) );
}

double
ProductionMechanism::generateMass(){
  
  if( m_type == kFlat ) return random( m_lowMass, m_highMass );
  
  double mass = 0;
  while( mass < m_lowMass || mass > m_highMass ){
    
    unsigned int channel = m_decGen();
    pair< double, double > bw = m_bwGen[channel]();
    
    mass = bw.first;
  }
  
  double prob = 0;
  for( unsigned int i = 0; i < m_bwGen.size(); ++i ){
    
    prob += m_bwGen[i].pdf( mass * mass ) * m_decGen.getProb( i );
  }
  
  // put in the factor of mass so resulting weights can be applied to 
  // obtain a distribution that is flat in mass instead of flat in s
  // (to get weights for reweighting flat in s remove the mass)
  
  m_lastWeight = 1 / ( prob * mass );
  
  return mass;
}

double
ProductionMechanism::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
ProductionMechanism::random( double low, double hi ) const {

        return( ( hi - low ) * gRandom->Uniform() + low );
}


