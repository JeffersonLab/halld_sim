/*
 *  GammaZToXYZ.cc
 *  GlueXTools
 *
 *  Modified GammaPToXYP.cc, replacing proton with heavy Z target
 *  Elton 4/14/2017
 *
 */

#include <iostream>
#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/GammaZToXYZ.h"
#include "AMPTOOLS_MCGEN/TwoBodyDecayFactory.h"

#include "IUAmpTools/Kinematics.h"

#include "UTILITIES/BeamProperties.h"

GammaZToXYZ::GammaZToXYZ( float lowMassXY, float highMassXY, 
                          float massX, float massY,
                          ProductionMechanism::Type type,
			  TString beamConfigFile  , float Bslope=2) : 

m_prodMech( ProductionMechanism::kZ, type, Bslope), // last arg is t dependence. Use value that is lower than any expected for reactions of interest. Elton 10/9/18
// m_target( 0, 0, 0, 108.),    // use mass of Tin
m_target( 0, 0, 0, 208.*0.931494),    // use mass of Pb since it is defined in particle tables.
m_childMass( 0 ) {

  m_childMass.push_back( massX );
  m_childMass.push_back( massY );
  
  m_prodMech.setMassRange( lowMassXY, highMassXY );
 
  // get beam properties from configuration file
  BeamProperties beamProp(beamConfigFile);
  cobrem_vs_E = (TH1D*)beamProp.GetFlux();
  cobrem_vs_E->GetName();
}

Kinematics* 
GammaZToXYZ::generate(){

  double beamE = cobrem_vs_E->GetRandom();
  m_beam.SetPxPyPzE(0,0,beamE,beamE);

  TLorentzVector resonance = m_prodMech.produceResonanceZ( m_beam);
  double genWeight = m_prodMech.getLastGeneratedWeight();
  
  vector< TLorentzVector > allPart;
  allPart.push_back( m_beam );
  // allPart.push_back( m_beam + m_target - resonance );
  
  TwoBodyDecayFactory decay( resonance.M(), m_childMass );
  
  vector<TLorentzVector> fsPart = decay.generateDecay();
  
  for( vector<TLorentzVector>::iterator aPart = fsPart.begin();
      aPart != fsPart.end(); ++aPart ){
    
    aPart->Boost( resonance.BoostVector() );
    allPart.push_back( *aPart );
  }

  allPart.push_back( m_beam + m_target - resonance );   // Move Recoil vector to position 3 after resonance
 
  return new Kinematics( allPart, genWeight );
}

void
GammaZToXYZ::addResonance( float mass, float width, float bf ){
  
  m_prodMech.addResonance( mass, width, bf );
}

