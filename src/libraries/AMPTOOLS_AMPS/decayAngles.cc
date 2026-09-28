//Goal: Calculate the decay angles for X->omegapi->4pi or Delta->protonpi decay, in either the helicity or Gottfried-Jackson reference frames. Inputs should be the relevant 4-momentum vectors in the lab frame, an integer flag specifying the desired reference frame, and a boolean specifying whether the decay happens at the upper or lower vertex of the t-channel reaction
//Gottfried-Jackson RF: The z-axis is equal to the direction of flight of the incoming beam photon in the parent rest frame.
//Adair RF: The z-axis is equal to the direction of flight of the incoming beam photon in the center of mass system.

#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>


#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "decayAngles.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

// Calculate the x, y, and z axes in the helicity reference frame
vector< TVector3 > getHelicityAxes(TVector3 parentCM, TVector3 inverseCM)
{
	TVector3 z = parentCM.Unit();
	TVector3 y = ( inverseCM.Unit() ).Cross( parentCM.Unit() ).Unit();
	TVector3 x = y.Cross( z );

	vector< TVector3 > xyz{x, y, z};
	return xyz;
} 

// Calculate the x, y, and z axes in the Gottfried-Jackson reference frame
vector< TVector3 > getGJAxes(TVector3 parentCM, TVector3 inverseCM, TVector3 inverseParent)
{
	TVector3 z = inverseParent.Unit();
	TVector3 y = ( inverseCM.Unit() ).Cross( parentCM.Unit() ).Unit();
	TVector3 x = y.Cross( z );

	vector< TVector3 > xyz{x, y, z};
	return xyz;
}

// Calculate production plane angle
double getPhiProd(double polAngle, TLorentzVector parentLab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex)
{
	// whichFrame = 1 for helicity, 2 for GJ
	assert( whichFrame == 1 || whichFrame == 2 );

	// Boost all P4 from lab to CM rest frame
	TLorentzVector cmLab = beamLab + targetLab;
	TVector3 boostCMToLab = cmLab.BoostVector();

	TLorentzVector parentCM = parentLab;
	TLorentzVector inverseLab;
	if( upperVertex == false )
		inverseLab = targetLab;
	else
		inverseLab = beamLab;
	TLorentzVector inverseCM = inverseLab;			
	parentCM.Boost( -1.0*boostCMToLab );
	inverseCM.Boost( -1.0*boostCMToLab );

	TVector3 boostParentToCM = parentCM.BoostVector();
	TLorentzVector inverseParent = inverseCM;
	inverseParent.Boost( -1.0*boostParentToCM );

	// convert P4s to P3
	TVector3 parentCM3 = parentCM.Vect();
	TVector3 inverseCM3 = inverseCM.Vect();
	TVector3 inverseParent3 = inverseParent.Vect();

	vector< TVector3 > locxyz;
	if( whichFrame == 1 )
		locxyz = getHelicityAxes(parentCM3, inverseCM3);
	else
		locxyz = getGJAxes(parentCM3, inverseCM3, inverseParent3);

	TVector3 y = locxyz[1];

	TVector3 eps( cos( polAngle ), sin( polAngle ), 0.0 );

	double phiProd = atan2( y.Dot( eps ), inverseCM.Vect().Unit().Dot( eps.Cross( y ) ) );

	return phiProd;	
} 


// Calculate angles for a one-step decay. If only one of the daughter particles has spin, that one should be used in this calculation 
vector< double > getOneStepAngles(TLorentzVector parentLab, TLorentzVector daughterLab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex)
{
	// whichFrame = 1 for helicity, 2 for GJ
	assert( whichFrame == 1 || whichFrame == 2 );

	// Boost all P4 from lab to CM rest frame
	TLorentzVector cmLab = beamLab + targetLab;
	TVector3 boostCMToLab = cmLab.BoostVector();

	TLorentzVector parentCM = parentLab;
	TLorentzVector daughterCM = daughterLab;
	TLorentzVector inverseCM;
	if( upperVertex == false )
		inverseCM = targetLab;
	else
		inverseCM = beamLab;
	parentCM.Boost( -1.0*boostCMToLab );
	daughterCM.Boost( -1.0*boostCMToLab );
	inverseCM.Boost( -1.0*boostCMToLab );

	// boost daughter and inverse to parent rest frame
	TVector3 boostParentToCM = parentCM.BoostVector();
	TLorentzVector daughterParent = daughterCM;
	TLorentzVector inverseParent = inverseCM;
	daughterParent.Boost( -1.0*boostParentToCM );
	inverseParent.Boost( -1.0*boostParentToCM );

	// convert P4s to P3
	TVector3 parentCM3 = parentCM.Vect();
	TVector3 inverseCM3 = inverseCM.Vect();
	TVector3 inverseParent3 = inverseParent.Vect();

	// define xyz axes in their own function (I think this will be useful later)
	vector< TVector3 > locxyz;
	if( whichFrame == 1 )
		locxyz = getHelicityAxes(parentCM3, inverseCM3);
	else
		locxyz = getGJAxes(parentCM3, inverseCM3, inverseParent3);

	TVector3 x = locxyz[0];
	TVector3 y = locxyz[1];
	TVector3 z = locxyz[2];

	// project daughter in parent's rest frame onto xyz axes
	TVector3 angles( daughterParent.Vect().Dot( x ),
				daughterParent.Vect().Dot( y ),
				daughterParent.Vect().Dot( z ) );

	double theta = angles.Theta();
	double phi = angles.Phi();

	vector< double > thetaPhi{theta, phi};

	return thetaPhi;
}


// Calculate angles for a two-step decay. This has been tested for photoproduction of a resonance decaying to omega and pi, with the omega subequently decaying to three pions. If the second step of the decay is a two-body decay, the input for granddaughter2Lab should be TLorentzVector(0,0,0,0)  
vector< double > getTwoStepAngles(TLorentzVector parentLab, TLorentzVector daughterLab, TLorentzVector granddaughter1Lab, TLorentzVector granddaughter2Lab, TLorentzVector beamLab, TLorentzVector targetLab, int whichFrame, bool upperVertex)
{
	// whichFrame = 1 for helicity, 2 for GJ
	assert( whichFrame == 1 || whichFrame == 2 );

	// Boost all P4 from lab to CM rest frame
	TLorentzVector cmLab = beamLab + targetLab;
	TVector3 boostCMToLab = cmLab.BoostVector();

	TLorentzVector parentCM = parentLab;
	TLorentzVector daughterCM = daughterLab;
	TLorentzVector granddaughter1CM = granddaughter1Lab;
	TLorentzVector granddaughter2CM = granddaughter2Lab;
	TLorentzVector inverseCM;
	if( upperVertex == false )
		inverseCM = targetLab;
	else
		inverseCM = beamLab;
	parentCM.Boost( -1.0*boostCMToLab );
	daughterCM.Boost( -1.0*boostCMToLab );
	granddaughter1CM.Boost( -1.0*boostCMToLab );
	granddaughter2CM.Boost( -1.0*boostCMToLab );
	inverseCM.Boost( -1.0*boostCMToLab );

	// boost daughter and inverse to parent rest frame
	TVector3 boostParentToCM = parentCM.BoostVector();
	TLorentzVector daughterParent = daughterCM;
	TLorentzVector granddaughter1Parent = granddaughter1CM;
	TLorentzVector granddaughter2Parent = granddaughter2CM;
	TLorentzVector inverseParent = inverseCM;
	daughterParent.Boost( -1.0*boostParentToCM );
	granddaughter1Parent.Boost( -1.0*boostParentToCM );
	granddaughter2Parent.Boost( -1.0*boostParentToCM );
	inverseParent.Boost( -1.0*boostParentToCM );

	// convert P4s to P3
	TVector3 parentCM3 = parentCM.Vect();
	TVector3 inverseCM3 = inverseCM.Vect();
	TVector3 inverseParent3 = inverseParent.Vect();

	// define xyz axes in their own function (I think this will be useful later)
	vector< TVector3 > locxyz;
	if( whichFrame == 1 )
		locxyz = getHelicityAxes(parentCM3, inverseCM3);
	else
		locxyz = getGJAxes(parentCM3, inverseCM3, inverseParent3);

	TVector3 x = locxyz[0];
	TVector3 y = locxyz[1];
	TVector3 z = locxyz[2];

	// project daughter in parent's rest frame onto xyz axes
	TVector3 daughterParent3 = daughterParent.Vect();
	TVector3 angles( daughterParent3.Dot( x ),
				daughterParent3.Dot( y ),
				daughterParent3.Dot( z ) );

	double theta = angles.Theta();
	double phi = angles.Phi();

	// Calculate decay of daughter particle in its helicity frame
	// Boost granddaughter(s) to daughter's rest frame
	TVector3 boostDaughterToParent = daughterParent.BoostVector();
	TLorentzVector granddaughter1Daughter = granddaughter1Parent;	
	TLorentzVector granddaughter2Daughter = granddaughter2Parent;
	granddaughter1Daughter.Boost( -1.0*boostDaughterToParent );	
	granddaughter2Daughter.Boost( -1.0*boostDaughterToParent );	


	vector< TVector3 > locxyzH = getHelicityAxes(daughterParent3, z);
	TVector3 xH = locxyzH[0];	
	TVector3 yH = locxyzH[1];	
	TVector3 zH = locxyzH[2];

	TVector3 daughterDecayVector;
	if( granddaughter2Lab.E() > 0) 
		daughterDecayVector = ( granddaughter1Daughter.Vect() ).Cross( granddaughter2Daughter.Vect() );
	else
		daughterDecayVector = granddaughter1Daughter.Vect();

	TVector3 anglesH( daughterDecayVector.Dot( xH ), 
				daughterDecayVector.Dot( yH ), 
				daughterDecayVector.Dot( zH ) );

	double thetaH = anglesH.Theta();
	double phiH = anglesH.Phi();	


  	// compute omega dalitz decay variable lambda
  	double m0 = 0.1349766;
  	double mq = 0.1395702;
  	double lambda_max = 3/4. * TMath::Power(1/9. * (5*daughterLab.M2() + 3*(m0*m0 - 4*mq*mq) - 4*sqrt(daughterLab.M2()*daughterLab.M2() + 3*daughterLab.M2()*(m0*m0-mq*mq))), 2); 
  	double lambda = fabs( daughterDecayVector.Dot( daughterDecayVector ) ) / lambda_max;

	vector< double > thetaPhiTwoStep{theta, phi, thetaH, phiH, lambda};

	return thetaPhiTwoStep;
}



