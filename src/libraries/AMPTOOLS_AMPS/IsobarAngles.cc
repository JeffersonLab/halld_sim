
#include <stdlib.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "AMPTOOLS_AMPS/IsobarAngles.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"

IsobarAngles::IsobarAngles( const vector< string >& args ) :
UserAmplitude< IsobarAngles >( args )
{

	assert( args.size() > 5 ); // at least 1 isobar
	
	m_jX      = atoi( args[0].c_str() ); // total J of produced resonance
	m_parX    = atoi( args[1].c_str() ); // parity of produced resonance
	m_lX      = atoi( args[2].c_str() ); // l between bachelor and isobar
	m_daughtX = string( args[3] );
	
	// loop over additional isobar parameters (J and daughters indices)
	m_nIsobars = 0;
	int maxPar = 4;
	for( unsigned int i = maxPar; i < args.size(); i+=2 ) {
		
		// total J of isobar
		m_jI.push_back( atoi( args[i].c_str() ) ); 
		
		// daughters of isobar (bachelor always first)
		m_daughtI.push_back( string( args[i+1] ) );  

		m_nIsobars++;
	}
	
	assert( m_jX >= 0  );
	assert( abs( (double)m_parX ) == 1 );
	assert( m_lX <= m_jX );
	for( int i = 0; i < m_jI[i]; i++) assert( m_jI[i] >= 0  );

}

complex< GDouble >
IsobarAngles::calcAmplitude( GDouble** pKin ) const
{

	TLorentzVector PX, Ptemp;
	TLorentzVector PIsobar[m_nIsobars], PBatch[m_nIsobars] ;
	pair<TLorentzVector, TLorentzVector> PNorm[m_nIsobars];
	
	for( unsigned int i = 0; i < m_daughtX.size(); ++i ){
		
		string num; num += m_daughtX[i];
		int index = atoi(num.c_str());
		Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
				  pKin[index][3], pKin[index][0] );
		PX += Ptemp;
	}
	
	for( int j = 0; j < m_nIsobars; j++ ){
		for( unsigned int i = 0; i < m_daughtI[j].size(); ++i ){
			
			string num; num += m_daughtI[j][i];
			int index = atoi(num.c_str());
			Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
					  pKin[index][3], pKin[index][0] );
			PIsobar[j] += Ptemp;
			if( i == 0 ) {
				PBatch[j] = Ptemp;
				PNorm[j].first = Ptemp; 
			}
			else if( i == 1 ) 
				PNorm[j].second = Ptemp;
			
		}
	}
	
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector PBatchX ( pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0] ); //= PX - PIsobar[0]; // batchelor from X decay
	
	TLorentzVector target ( 0, 0, 0, 0.938 );
	TLorentzVector cms = beam + target;

	//
	TVector3 cmsBoost = cms.BoostVector();
	TLorentzVector recoilCms = recoil;
	recoilCms.Boost(-1.0*cmsBoost);
	
	// calculate decay angles in resonance X rest frame
	TVector3 XRestBoost = PX.BoostVector();
	
	TLorentzVector beamX   = beam;
	TLorentzVector recoilX = recoil;
	TLorentzVector batchX  = PBatchX;
	beamX.Boost(-1.0*XRestBoost);
	recoilX.Boost(-1.0*XRestBoost);
	batchX.Boost(-1.0*XRestBoost);
	
	TVector3 z = -recoilX.Vect().Unit();
	TVector3 y = (beam.Vect().Unit()).Cross(z).Unit();
	TVector3 x = y.Cross(z);
	//

	/*
	TLorentzRotation cmsBoost( -cms.BoostVector() );
	TLorentzVector recoilCms = cmsBoost * recoil;

	// calculate decay angles in resonance X rest frame
	TLorentzRotation XRestBoost( -PX.BoostVector() );
	
	TLorentzVector beamX   = XRestBoost * beam;
	TLorentzVector batchX  = XRestBoost * PBatchX;
	
	TVector3 z = -recoilCms.Vect().Unit();
	TVector3 y = beamX.Vect().Cross(z).Unit();
	TVector3 x = y.Cross(z);
	*/	

	TVector3 anglesBatchX( (batchX.Vect()).Dot(x),
			       (batchX.Vect()).Dot(y),
			       (batchX.Vect()).Dot(z) );
	
	GDouble cosThetaBatchX = anglesBatchX.CosTheta();
	GDouble phiBatchX = anglesBatchX.Phi();
	
	// calculate decay angles in isobar rest frame
	vector<GDouble> cosThetaIso, phiIso, k, q;
	TVector3 zIsoPrevious = z;
	for( int i = 0; i < m_nIsobars; i++ ){
		
		TLorentzRotation isoRestBoost( -PIsobar[i].BoostVector() );
		TLorentzVector PBatchIso = isoRestBoost * PBatch[i];
		TLorentzVector PResonanceIso = isoRestBoost * (PIsobar[i] - PBatch[i]);
		TLorentzVector PNormIso1 = isoRestBoost * PNorm[i].first;
		TLorentzVector PNormIso2 = isoRestBoost * PNorm[i].second;

		TVector3 zIso = PResonanceIso.Vect().Unit();
		// only true for b1 pi otherwise need original z for omega-rho...
		TVector3 yIso = zIsoPrevious.Cross(zIso).Unit(); 
		TVector3 xIso = yIso.Cross(zIso);
		
		TVector3 PAngles = PBatchIso.Vect();
		if(m_daughtI[i].size() == 3 and m_daughtI.size() == uint(m_nIsobars)) // 3-body decay (e.g. omega) 
			PAngles = (PNormIso1.Vect()).Cross(PNormIso2.Vect());

		// why uses axes in X rest frame?  Is this GJ or helicity?
		TVector3 anglesIso( (PAngles).Dot(xIso),
				    (PAngles).Dot(yIso),
				    (PAngles).Dot(zIso) );
		
		cosThetaIso.push_back(anglesIso.CosTheta());
		phiIso.push_back(anglesIso.Phi());
		
		k.push_back( breakupMomentum( PX.M(), PIsobar[i].M(), PBatchX.M() ) );
		q.push_back( breakupMomentum( PIsobar[i].M(), PBatch[i].M(), (PIsobar[i] - PBatch[i]).M() ) );
		
		zIsoPrevious = zIso;
	}  
	
	const vector< int >& perm = getCurrentPermutation();
	
	complex< GDouble > i( 0, 1 );
	complex< GDouble > ans( 0, 0 );
	
	// in general we also need a sum over resonance helicities here
	// however, we assume a production mechanism that only produces
	// resonance helicities +-1
	
	for( int mL = -m_lX; mL <= m_lX; ++mL ){
		
		complex< GDouble > term( 1, 0 );

		/*
		for ( int i = 0; i < m_nIsobars; i++ ){
			
			for( int mI = -m_jI[i]; mI <= m_jI[i]; ++mI ){
				
				//fixed helicity = +1 for now
				term += Y( m_jI[i], mI, cosThetaIso[i], phiIso[i] );
				cout<<"Isobar "<<i<<": "<<m_jI[i]<<" "<<mI<<" "<<cosThetaIso[i]<<" "<<phiIso[i]<<" "<<Y( m_jI[i], mI, cosThetaIso[i], phiIso[i] )<<endl;
				// *clebschGordan( m_jI[i], m_lX, mI, mL, m_jX, 1 ); 
			}	
		}
		*/

		term *= Y( m_lX, mL, cosThetaBatchX, phiBatchX );
		//cout<<"X: "<<m_lX<<" "<<mL<<" "<<cosThetaBatchX<<" "<<phiBatchX<<" "<<Y( m_lX, mL, cosThetaBatchX, phiBatchX )<<endl;
		ans += term;
	}
	//cout<<"Total amplitude = "<<ans<<endl;

	ans *= cos(recoil.Phi());
	ans *= cos(phiBatchX);
	
	return ans;
}
