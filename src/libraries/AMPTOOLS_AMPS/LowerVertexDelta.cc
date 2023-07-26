
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

//	cout << "cos(theta_p)_GJ = " << userVars[kCosTheta] << ", phi_p_GJ = " << userVars[kPhi] << endl;
/*
	// boost all 4-vectors to the COM (gammap) frame
	TLorentzVector gammap = target + beam;
	TLorentzRotation gammapBoost( -gammap.BoostVector() );
	TLorentzVector beam_gammapRF = gammapBoost * beam;
	TLorentzVector target_gammapRF = gammapBoost * target;
	TLorentzVector p1_gammapRF = gammapBoost * p1;
	TLorentzVector p2_gammapRF = gammapBoost * p2;
	TLorentzVector p3_gammapRF = gammapBoost * p3;
	

	TLorentzVector lowerVertexResonance = p1_gammapRF + p2_gammapRF;
	TLorentzRotation lowerVertexBoost( -lowerVertexResonance.BoostVector() );

	TLorentzVector target_lowerVertexRF = lowerVertexBoost * target_gammapRF;
	TLorentzVector beam_lowerVertexRF = lowerVertexBoost * beam_gammapRF;
	TLorentzVector upperVertexResonance_lowerVertexRF = lowerVertexBoost * p3_gammapRF;
	TLorentzVector p2_lowerVertexRF = lowerVertexBoost * p2_gammapRF;
	TLorentzVector p1_lowerVertexRF = lowerVertexBoost * p1_gammapRF;
	
	// normal to the production plane
	TVector3 y = (target_lowerVertexRF.Vect().Unit().Cross(upperVertexResonance_lowerVertexRF.Vect().Unit())).Unit();
	// choose Gottfried-Jackson frame: z-axis along -target direction in baryon rest frame
	TVector3 z = target_lowerVertexRF.Vect().Unit();
	TVector3 x = y.Cross(z).Unit();

	// test with GJ angles:
	TVector3 z = target_lowerVertexRF.Vect().Unit();
	TVector3 y = (upperVertexResonance_lowerVertexRF.Vect()).Cross(z).Unit();
	TVector3 x = y.Cross(z).Unit();
	
	TVector3 angles( (p2_lowerVertexRF.Vect()).Dot(x),
					(p2_lowerVertexRF.Vect()).Dot(y),
					(p2_lowerVertexRF.Vect()).Dot(z) );

	userVars[kCosTheta]	= angles.CosTheta();
	userVars[kPhi] 		= angles.Phi();

	TLorentzVector pDelta = p1 + p2;
	TLorentzVector ptot_DeltaRF = lowerVertexBoost * lowerVertexResonance;

	cout << "cos(theta_p) = " << userVars[kCosTheta] << ", phi_p = " << userVars[kPhi] << endl;
	cout << "In the lab frame: " << endl;
	cout << "Delta: ";
	pDelta.Print();
	cout << "Target: ";
	target.Print();
	
	cout << "In CM (gammap) frame: " << endl;
	cout << "Delta: ";
	lowerVertexResonance.Print();
	cout << "Target: ";
	target_gammapRF.Print();

	cout << "In Delta rest frame: " << endl;
	cout << "Delta: ";
	ptot_DeltaRF.Print();
	cout << "Target: ";
	target_lowerVertexRF.Print();
	cout << "Beam: ";
	beam_lowerVertexRF.Print();
	cout << "omegapi: ";
	upperVertexResonance_lowerVertexRF.Print();
	cout << "Recoil pion: ";
	p2_lowerVertexRF.Print();
	cout << "Recoil proton: ";
	p1_lowerVertexRF.Print();
	

	cout << "GJ z-axis: "; 
	z.Print();
	cout << "GJ y-axis: "; 
	y.Print();
	cout << "GJ x-axis: "; 
	x.Print();
*/
}

#ifdef GPU_ACCELERATION
void
LowerVertexDelta::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPULowerVertexDelta_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
			m_d, m_p, m_c, m_s );
}

#endif // GPU_ACCELERATION
