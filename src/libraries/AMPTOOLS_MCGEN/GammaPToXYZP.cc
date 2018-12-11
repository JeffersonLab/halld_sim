/*
 *  GammaPToXYZP.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 5/25/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "GammaPToXYZP.h"
#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/DalitzDecayFactory.h"
#include "IUAmpTools/Kinematics.h"

#include "UTILITIES/BeamProperties.h"

GammaPToXYZP::GammaPToXYZP( float lowMassXYZ, float highMassXYZ, 
        float massX, float massY, float massZ,
        ProductionMechanism::Type type,
        float slope,
	TString beamConfigFile ) :
    m_prodMech( ProductionMechanism::kProton, type, slope), // last arg is t dependence
    m_target( 0, 0, 0, 0.938 ),
    m_childMass( 0 ) 
{

    m_childMass.push_back( massX );
    m_childMass.push_back( massY );
    m_childMass.push_back( massZ );

    m_prodMech.setMassRange( lowMassXYZ, highMassXYZ );

    // get beam properties from configuration file
    BeamProperties beamProp(beamConfigFile);
    cobrem_vs_E = (TH1D*)beamProp.GetFlux();
    cobrem_vs_E->GetName();
}

Kinematics* 
GammaPToXYZP::generate(){

    double beamE = cobrem_vs_E->GetRandom();
    m_beam.SetPxPyPzE(0,0,beamE,beamE);

    TLorentzVector resonance = m_prodMech.produceResonance( m_beam );
    double genWeight = m_prodMech.getLastGeneratedWeight();

    vector< TLorentzVector > allPart;
    allPart.push_back( m_beam );
    allPart.push_back( m_beam + m_target - resonance );

    DalitzDecayFactory decay( resonance.M(), m_childMass );

    vector<TLorentzVector> fsPart = decay.generateDecay();

    for( vector<TLorentzVector>::iterator aPart = fsPart.begin();
            aPart != fsPart.end(); ++aPart ){

        aPart->Boost( resonance.BoostVector() );
        allPart.push_back( *aPart );
    }

    return new Kinematics( allPart, genWeight );
}

void
GammaPToXYZP::addResonance( float mass, float width, float bf ){

    m_prodMech.addResonance( mass, width, bf );
}

