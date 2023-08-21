
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/LowerVertexDelta.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/decayAngles.h"

LowerVertexDelta::LowerVertexDelta( const vector< string >& args ) :
UserAmplitude< LowerVertexDelta >( args )
{
	assert( args.size() == 5 );

	m_d = atoi( args[0].c_str() ); // Twice the helicity of decaying Delta baryon: (3,1,-1, or -3)
	m_p = atoi( args[1].c_str() ); // Twice the helicity of final state proton (+/-1)
	m_c = atoi( args[2].c_str() ); // Wigner D function (+1) or its complex conjugate (-1)
	m_s = atoi( args[3].c_str() ); // The amplitude gets multiplied by a positive or negative sign

	lowerVertex = args[4].c_str(); // indices of proton and pi+ from the lower vertex

	// make sure values are reasonable
	assert( abs( m_d ) == 1 || abs( m_d ) == 3 );
	assert( abs( m_p ) == 1 );
	assert( abs( m_c ) == 1 );
	assert( abs( m_s ) == 1 );
}


complex< GDouble >
LowerVertexDelta::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {	

	GDouble phi		= userVars[kPhi];
	GDouble cosTheta	= userVars[kCosTheta];

	GDouble lambda_Delta 	= m_d / 2.;
	GDouble lambda_proton 	= m_p / 2.;

	complex <GDouble> amplitude;
	if( m_c == 1 ){
		amplitude = 2/3. * m_s * wignerD( 3/2., lambda_Delta, lambda_proton, cosTheta, phi );
	}
	else{
		amplitude = 2/3. * m_s * conj( wignerD( 3/2., lambda_Delta, lambda_proton, cosTheta, phi ) );
	}

	return amplitude;		
}

void
LowerVertexDelta::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

	TLorentzVector target ( 0, 0, 0, 0.9382720813);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector p1, p2, pDelta;
	
	string lv1; lv1 += lowerVertex[0];
	string lv2; lv2 += lowerVertex[1];

        int index1 = atoi( lv1.c_str() );     
	int index2 = atoi( lv2.c_str() );

	p1.SetPxPyPzE( pKin[index1][1], pKin[index1][2], pKin[index1][3], pKin[index1][0] );
	p2.SetPxPyPzE( pKin[index2][1], pKin[index2][2], pKin[index2][3], pKin[index2][0] );

	pDelta = p1 + p2;
	
	vector< double > thetaPhi = getOneStepAngles( pDelta, p1, beam, target, 2, false );

	userVars[kCosTheta] 	= TMath::Cos( thetaPhi[0] );
	userVars[kPhi]		= thetaPhi[1];

	return;
}

#ifdef GPU_ACCELERATION
void
LowerVertexDelta::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPULowerVertexDelta_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
			m_d, m_p, m_c, m_s );
}

#endif // GPU_ACCELERATION
