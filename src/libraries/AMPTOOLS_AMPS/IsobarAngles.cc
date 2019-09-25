
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

		m_nIsobars++;
	}
	
	assert( m_jX >= 0  );
	assert( m_lX <= m_jX );
	for( unsigned int i = 0; i < m_jI.size(); i++) {
		cout<<"Isobar: J="<<m_jI[i]<<" and daughters="<<m_daughtI[i].data()<<endl;
		assert( m_jI[i] >= 0  );
		
		// check if this is a subseqent bachelor decay
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
	//GDouble pTest[7][4] = {{8.605265,0.000000,0.000000,8.605265}, {0.957474,-0.016470,-0.008349,0.189899}, {1.698119,-0.188228,0.001651,1.681872}, {1.903256,-0.357289,0.318952,1.836714}, {0.797837,0.207803,0.013740,0.757425}, {2.544201,0.121276,-0.092830,2.535775}, {1.642650,0.232908,-0.233165,1.603580}};

	TLorentzVector PX, Ptemp;
	TLorentzVector PIsobar[m_nIsobars], PBach[m_nIsobars];
	TLorentzVector PIsobarX[m_nIsobars], PBachX[m_nIsobars];
	pair<TLorentzVector, TLorentzVector> PNorm[m_nIsobars], PNormX[m_nIsobars];
	
	// add particle P4s to get momentum of X
	for( unsigned int i = 0; i < m_daughtX.size(); ++i ){
		
		string num; num += m_daughtX[i];
		int index = atoi(num.c_str());
		Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
				  pKin[index][3], pKin[index][0] );
		//Ptemp.SetPxPyPzE( pTest[index][1], pTest[index][2],
		//		  pTest[index][3], pTest[index][0] );
		PX += Ptemp;
	}
	
	// add particle P4s to get momentum of Isobars
	for( int j = 0; j < m_nIsobars; j++ ){
		for( unsigned int i = 0; i < m_daughtI[j].size(); ++i ){
			
			string num; num += m_daughtI[j][i];
			int index = atoi(num.c_str());
			Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
					  pKin[index][3], pKin[index][0] );
			//Ptemp.SetPxPyPzE( pTest[index][1], pTest[index][2],
			//		  pTest[index][3], pTest[index][0] );
			PIsobar[j] += Ptemp;
			if( i == 0 ) {
				PBach[j] = Ptemp;
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
	TLorentzVector PBachXdecay = PX - PIsobar[0]; // bachelor from X decay
	
	// calculate decay angles in resonance X rest frame
	TVector3 XRestBoost = PX.BoostVector();
	
	TLorentzVector beamX   = beam;
	TLorentzVector recoilX = recoil;
	TLorentzVector isobarX = PIsobar[0];
	beamX.Boost(-1.0*XRestBoost);
	recoilX.Boost(-1.0*XRestBoost);
	isobarX.Boost(-1.0*XRestBoost);

	// keep vectors for isobars in X rest frame for later angle definitions
	for( int i = 0; i < m_nIsobars; i++ ){
		TLorentzVector temp;
		temp = PIsobar[i]; temp.Boost(-1.0*XRestBoost); PIsobarX[i] = temp;
		temp = PBach[i]; temp.Boost(-1.0*XRestBoost); PBachX[i] = temp;
		temp = PNorm[i].first; temp.Boost(-1.0*XRestBoost); PNormX[i].first = temp;
		temp = PNorm[i].second; temp.Boost(-1.0*XRestBoost); PNormX[i].second = temp;
	}

	// For GJ frame: choose beam as z-axis for reference 
	TVector3 z = beamX.Vect().Unit();
	TVector3 y = (beamX.Vect().Unit()).Cross((-recoilX.Vect().Unit())).Unit();
	TVector3 x = y.Cross(z);
	
	TVector3 anglesX( (isobarX.Vect()).Dot(x),
			  (isobarX.Vect()).Dot(y),
			  (isobarX.Vect()).Dot(z) );
	
	GDouble cosThetaX = anglesX.CosTheta();
	GDouble phiX = anglesX.Phi();
	
	////////////////////////////////////////////////////////////
	// calculate decay angles in isobar rest (helicity) frame //
	////////////////////////////////////////////////////////////
	vector<GDouble> cosThetaIso, phiIso, k, q;
	pair<TLorentzVector, TLorentzVector> isoPrevious;
	for( int i = 0; i < m_nIsobars; i++ ){
		
		// boost from X rest frame to isobar rest frame (could do in terms of previous frame)
		TVector3 isoRestBoost = PIsobarX[i].BoostVector();
		TLorentzVector PResonanceIso = PIsobarX[i] - PBachX[i];
		TLorentzVector PNormIso1 = PNormX[i].first;
		TLorentzVector PNormIso2 = PNormX[i].second;
		PResonanceIso.Boost(-1.0*isoRestBoost);
		PNormIso1.Boost(-1.0*isoRestBoost);
		PNormIso2.Boost(-1.0*isoRestBoost);

		// Helicity frame z-axis is direction of isobar in X rest frame by default
		TVector3 zIso = PIsobarX[i].Vect().Unit(); 
		TVector3 yIso = (z.Cross(zIso)).Unit(); // decay plane from X rest frame
		
		// later stage of single bachelor decays (eg. omega->3pi in b1pi or K*->Kpi in K1K production)
		if(m_isBach[i]) { 

			// boost from X frame to previous isobar frame
			TVector3 boost1 = isoPrevious.first.BoostVector();
			// boost to previous isobar frame to current isobar frame
			TVector3 boost2 = isoPrevious.second.BoostVector(); 
			
			PResonanceIso = PIsobarX[i] - PBachX[i];
			PResonanceIso.Boost(-1.0*boost1); PResonanceIso.Boost(-1.0*boost2);
			PNormIso1 = PNormX[i].first; PNormIso2 = PNormX[i].second;
			PNormIso1.Boost(-1.0*boost1); PNormIso2.Boost(-1.0*boost1); 
			PNormIso1.Boost(-1.0*boost2); PNormIso2.Boost(-1.0*boost2);
			
			zIso = (isoPrevious.second.Vect()).Unit();
			yIso = ((isoPrevious.first.Vect().Unit()).Cross(zIso)).Unit();
		}
		TVector3 xIso = yIso.Cross(zIso);
		
		TVector3 PAngles = PResonanceIso.Vect();
		if(m_daughtI[i].size() == 3 && i+1 == m_nIsobars) // 3-body decay (e.g. omega) 
			PAngles = ((PNormIso1.Vect()).Cross(PNormIso2.Vect())); 
		
		// Angles in isobar rest frame
		TVector3 anglesIso( (PAngles).Dot(xIso),
				    (PAngles).Dot(yIso),
				    (PAngles).Dot(zIso) );
		
		cosThetaIso.push_back(anglesIso.CosTheta());
		phiIso.push_back(anglesIso.Phi());

		k.push_back( breakupMomentum( PX.M(), PIsobar[i].M(), PBachXdecay.M() ) );
		q.push_back( breakupMomentum( PIsobar[i].M(), PBach[i].M(), (PIsobar[i] - PBach[i]).M() ) );
		
		// reference vector for later step in previous frame
		isoPrevious.first = PIsobarX[i];
		// reference vector for later step in current frame
		isoPrevious.second = PResonanceIso; 
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
				
				//cout<<"Isobar "<<i<<": "<<m_jI[i]<<" "<<mI<<" "<<cosThetaIso[i]<<" "<<phiIso[i]<<" "<<Y( m_jI[i], mI, cosThetaIso[i], phiIso[i] )<<" "<<sin(acos(cosThetaIso[i]))<<endl;
				
				// fixed helicity = +1 in CG (to be multiplied to termIso later)
				//* clebschGordan( m_jI[i], m_lX, mI, mL, m_jX, 1 );
			}	
			
			if( i==0 ) term = termIso; // add 1st isobar sum over mI to the given mL term
			else term *= termIso; // multiply additional isobar sum to the same mL term
		}

		// decay angle constribution for X decay
		term *= Y( m_lX, mL, cosThetaX, phiX );

		// add each mL term to total amplitude
		ans += term;
	}
	
	return ans;
}
