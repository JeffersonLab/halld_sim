
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

	assert( args.size() > 4 ); // at least 1 isobar
	
	m_jX      = atoi( args[0].c_str() ); // total J of produced resonance
	m_lX      = atoi( args[1].c_str() ); // l between bachelor and isobar
	m_daughtX = string( args[2] );
	
	// loop over additional isobar parameters (J and daughters indices)
	m_nIsobars = 0;
	int maxPar = 3;
	for( unsigned int i = maxPar; i < args.size(); i+=2 ) {
		
		// total J of isobar
		m_jI.push_back( atoi( args[i].c_str() ) ); 
		
		// daughters of isobar (bachelor always first)
		m_daughtI.push_back( string( args[i+1] ) );  

		/*
		// check if this is a subseqent batchelor decay
		bool isBach = false;
		if(m_daughtI.size() > 1) {
			if(m_daughtI[i-maxPar].size() == (m_daughtI[i-maxPar-1].size() - 1)) {
				isBach = true;
				cout<<"Found subsequent bachelor decay "<<m_daughtI[i-maxPar-1].data()<<" -> "<<m_daughtI[i-maxPar].data()<<endl;
			}
		}
		m_isBach.push_back( isBach );
		*/

		m_nIsobars++;
	}
	
	assert( m_jX >= 0  );
	assert( m_lX <= m_jX );
	for( unsigned int i = 0; i < m_jI.size(); i++) {
		cout<<"Isobar: J="<<m_jI[i]<<" and daughters="<<m_daughtI[i].data()<<endl;
		assert( m_jI[i] >= 0  );
		
		// check if this is a subseqent batchelor decay
		bool isBach = false;
		if(i>0 && m_daughtI[i].size() == (m_daughtI[i-1].size() - 1)) {
			isBach = true;
			cout<<"Found subsequent bachelor decay "<<m_daughtI[i-1].data()<<" -> "<<m_daughtI[i].data()<<endl;
		}
		m_isBach.push_back( isBach );
	}
}

complex< GDouble >
IsobarAngles::calcAmplitude( GDouble** pKin ) const
{

	TLorentzVector PX, Ptemp;
	TLorentzVector PIsobar[m_nIsobars], PBatch[m_nIsobars];
	TLorentzVector PIsobarX[m_nIsobars], PBatchX[m_nIsobars];
	pair<TLorentzVector, TLorentzVector> PNorm[m_nIsobars], PNormX[m_nIsobars];
	
	// add particle P4s to get momentum of X
	for( unsigned int i = 0; i < m_daughtX.size(); ++i ){
		
		string num; num += m_daughtX[i];
		int index = atoi(num.c_str());
		Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
				  pKin[index][3], pKin[index][0] );
		PX += Ptemp;
	}
	
	// add particle P4s to get momentum of Isobars
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
	
	/////////////////////////////////
	// calculate decay angles of X //
	/////////////////////////////////
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector PBatchXdecay = PX - PIsobar[0]; // batchelor from X decay
	
	// calculate decay angles in resonance X rest frame
	TVector3 XRestBoost = PX.BoostVector();
	
	TLorentzVector beamX   = beam;
	TLorentzVector recoilX = recoil;
	TLorentzVector batchX  = PBatchXdecay;
	beamX.Boost(-1.0*XRestBoost);
	recoilX.Boost(-1.0*XRestBoost);
	batchX.Boost(-1.0*XRestBoost);

	// keep vectors for isobars in X rest frame for later angle definitions
	for( int i = 0; i < m_nIsobars; i++ ){
		TLorentzVector temp;
		temp = PIsobar[i]; temp.Boost(-1.0*XRestBoost); PIsobarX[i] = temp;
		temp = PBatch[i]; temp.Boost(-1.0*XRestBoost); PBatchX[i] = temp;
		temp = PNormX[i].first; temp.Boost(-1.0*XRestBoost); PNormX[i].first = temp;
		temp = PNormX[i].second; temp.Boost(-1.0*XRestBoost); PNormX[i].second = temp;
	}
	
	TVector3 z = beamX.Vect().Unit();
	TVector3 y = (beamX.Vect().Unit()).Cross((-recoilX.Vect().Unit())).Unit();
	TVector3 x = y.Cross(z);
	
	TVector3 anglesBatchX( (batchX.Vect()).Dot(x),
			       (batchX.Vect()).Dot(y),
			       (batchX.Vect()).Dot(z) );
	
	GDouble cosThetaBatchX = anglesBatchX.CosTheta();
	GDouble phiBatchX = anglesBatchX.Phi();
	
	/////////////////////////////////////////////////
	// calculate decay angles in isobar rest frame //
	/////////////////////////////////////////////////
	vector<GDouble> cosThetaIso, phiIso, k, q;
	pair<TVector3, TVector3> zIsoPrevious;
	for( int i = 0; i < m_nIsobars; i++ ){
		
		// boost from X rest frame to isobar rest frame (could do in terms of previous frame)
		TVector3 isoRestBoost = PIsobarX[i].BoostVector();
		TLorentzVector PBatchIso = PBatchX[i];
		TLorentzVector PResonanceIso = PIsobarX[i] - PBatchX[i];
		TLorentzVector PNormIso1 = PNormX[i].first;
		TLorentzVector PNormIso2 = PNormX[i].second;
		PBatchIso.Boost(-1.0*isoRestBoost);
		PResonanceIso.Boost(-1.0*isoRestBoost);
		PNormIso1.Boost(-1.0*isoRestBoost);
		PNormIso2.Boost(-1.0*isoRestBoost);

		// Helicity frame z-axis is direction of isobar in X rest frame by default
		TVector3 zIso = PIsobarX[i].Vect().Unit(); 
		TVector3 yIso = (z.Cross(zIso)).Unit(); // decay plane from X rest frame
		
		// later stage of single batchelor decays (eg. omega->3pi in b1pi production)
		if(m_isBach[i]) { 
			zIso = zIsoPrevious.first;
			yIso = zIsoPrevious.second.Cross(zIsoPrevious.first); 
		}
		TVector3 xIso = yIso.Cross(zIso);
		
		TVector3 PAngles = PBatchIso.Vect();
		if(m_daughtI[i].size() == 3 and m_daughtI.size() == uint(m_nIsobars)) // 3-body decay (e.g. omega) 
			PAngles = (PNormIso1.Vect()).Cross(PNormIso2.Vect());
		
		// Angles in isobar rest frame
		TVector3 anglesIso( (PAngles).Dot(xIso),
				    (PAngles).Dot(yIso),
				    (PAngles).Dot(zIso) );
		
		cosThetaIso.push_back(anglesIso.CosTheta());
		phiIso.push_back(anglesIso.Phi());

		k.push_back( breakupMomentum( PX.M(), PIsobar[i].M(), PBatchXdecay.M() ) );
		q.push_back( breakupMomentum( PIsobar[i].M(), PBatch[i].M(), (PIsobar[i] - PBatch[i]).M() ) );
		
		/////////// Need to keep P4 here instead of z-axis ///////////////
		// reference vector for later step in current frame
		zIsoPrevious.first = PAngles.Unit(); 
		// reference vector for later step in previous frame
		zIsoPrevious.second = zIso;
	}  
	
	const vector< int >& perm = getCurrentPermutation();
	
	complex< GDouble > i( 0, 1 );
	complex< GDouble > ans( 0, 0 ); // total amplitude to be returned 
		
	// loop over possible orbital angular momentum (mL) in X decay
	for( int mL = -m_lX; mL <= m_lX; ++mL ){
		
		// can fix mL by hand to explore different Ylm's
		if(mL != 0) continue;

		complex< GDouble > term( 0, 0 );

		// loop over possible isobars given in config file and calculate contribution to amplitude
		for ( int i = 0; i < m_nIsobars; i++ ){
			
			complex< GDouble > termIso( 0, 0 );
			
			// loop over possible orbital angular momentum (mI) in Isobar decay
			for( int mI = -m_jI[i]; mI <= m_jI[i]; ++mI ){
	
				if(mI != 0) continue;

				// decay angle contribution for Isobar decay (replace with wignerD later)
				termIso += Y( m_jI[i], mI, cosThetaIso[i], phiIso[i] ); 
				
				//cout<<"Isobar "<<i<<": "<<m_jI[i]<<" "<<mI<<" "<<cosThetaIso[i]<<" "<<phiIso[i]<<" "<<Y( m_jI[i], mI, cosThetaIso[i], phiIso[i] )<<endl;
				
				// fixed helicity = +1 in CG (to be multiplied to termIso later)
				//* clebschGordan( m_jI[i], m_lX, mI, mL, m_jX, 1 );
			}	
			
			if( i==0 ) term = termIso; // add 1st isobar sum over mI to the given mL term
			else term *= termIso; // multiply additional isobar sum to the same mL term
		}

		// decay angle constribution for X decay
		term *= Y( m_lX, mL, cosThetaBatchX, phiBatchX );

		// add each mL term to total amplitude
		ans += term;
	}
	
	return ans;
}
