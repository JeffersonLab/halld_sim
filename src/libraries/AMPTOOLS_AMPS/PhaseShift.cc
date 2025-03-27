#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/PhaseShift.h"

PhaseShift::PhaseShift( const vector< string >& args ) :
UserAmplitude< PhaseShift >( args )
{
    assert( args.size() == 4 );

    m_daughters = pair< string, string >( args[0], args[1] );

    m_p0 = AmpParameter( args[2] );
    m_p1 = AmpParameter( args[3] );

    registerParameter( m_p0 );
    registerParameter( m_p1 );
}

complex< GDouble >
PhaseShift::calcAmplitude( GDouble** pKin ) const
{
    TLorentzVector p1, p2, ptot, ptemp;

    for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
        string num; num += m_daughters.first[i];
        int index = atoi( num.c_str() );
        ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0] );

        p1 += ptemp;
        ptot += ptemp;
    }

    for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
        string num; num += m_daughters.second[i];
        int index = atoi( num.c_str() );
        ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0] );

        p2 += ptemp;
        ptot += ptemp;
    }

    double q = fabs( breakupMomentum( ptot.M(), p1.M(), p2.M() ) );

    double deltaJ = m_p0 + m_p1*q; // some function of q

    complex< GDouble > fJ = ( TMath::Cos( deltaJ ) , TMath::Sin( deltaJ ) );

    fJ *= TMath::Sin( deltaJ ) / q;

    return fJ;

}
