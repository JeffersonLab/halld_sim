
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/DeltaAngles_res.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

DeltaAngles_res::DeltaAngles_res( const vector< string >& args ) :
UserAmplitude< DeltaAngles_res >( args )
{
	assert( args.size() == 13 || args.size() == 12 );
	
	rho011  = AmpParameter( args[0] );
	rho031  = AmpParameter( args[1] );
	rho03m1 = AmpParameter( args[2] );
	
	rho111  = AmpParameter( args[3] );
	rho133  = AmpParameter( args[4] );
	rho131  = AmpParameter( args[5] );
	rho13m1 = AmpParameter( args[6] );
	
	rho231  = AmpParameter( args[7] );
	rho23m1 = AmpParameter( args[8] );

	lowerVertex = args[9].c_str();
	upperVertex = args[10].c_str();
	
	// need to register any free parameters so the framework knows about them
	registerParameter( rho011 );
	registerParameter( rho031 );
	registerParameter( rho03m1 );
	
	registerParameter( rho111 );
	registerParameter( rho133 );
	registerParameter( rho131 );
	registerParameter( rho13m1 );
	
	registerParameter( rho231 );
	registerParameter( rho23m1 );

	if(args.size() == 13){
		polAngle  = atof(args[11].c_str() ); // azimuthal angle of the photon polarization vector in the lab.
		polFraction = AmpParameter( args[12] ); // fraction of polarization (0-1)
		std::cout << "Fixed polarisation of " << polFraction << " and angle of " << polAngle << " degrees." << std::endl;
	}
	else
		assert(0);
}


complex< GDouble >
DeltaAngles_res::calcAmplitude( GDouble** pKin ) const {
	
	TLorentzVector target ( 0, 0, 0, 0.9382720813);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector p1, p2, p3, ptot, ptemp; //p1 = proton from decaying lower vertex, p2 = pi+ from decaying lower vertex, p3 = upper vertex resonance (b1 in this case)

	//lower vertex proton is always index 1, pi+ is always index 2

	p1.SetPxPyPzE( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	p2.SetPxPyPzE( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );

	ptot = p1 + p2;

	//p3 is sum of all particles in upper vertex
	for( unsigned int i = 0; i < upperVertex.size(); ++i ){
		string num; num += upperVertex[i];
		int index = atoi(num.c_str());
		ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0] );
		p3 += ptemp;
		ptot += ptemp;
	}

	TLorentzVector lowerVertexResonance = p1 + p2;
	TLorentzRotation lowerVertexBoost( -lowerVertexResonance.BoostVector() );

	TLorentzVector target_lowerVertexRF = lowerVertexBoost * target;
	TLorentzVector beam_lowerVertexRF = lowerVertexBoost * beam;
	TLorentzVector upperVertexResonance_lowerVertexRF = lowerVertexBoost * p3;
	TLorentzVector p2_lowerVertexRF = lowerVertexBoost * p2;
	
	// normal to the production plane
	TVector3 y = (target_lowerVertexRF.Vect().Unit().Cross(upperVertexResonance_lowerVertexRF.Vect().Unit())).Unit();
	// choose Gottfried-Jackson frame: z-axis along -target direction in baryon rest frame
	TVector3 z = target_lowerVertexRF.Vect().Unit();
	TVector3 x = y.Cross(z).Unit();
	
	TVector3 angles( (p2_lowerVertexRF.Vect()).Dot(x),
					(p2_lowerVertexRF.Vect()).Dot(y),
					(p2_lowerVertexRF.Vect()).Dot(z) );
	
	GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
	GDouble cosSqTheta = cos(angles.Theta())*cos(angles.Theta());
	GDouble sin2Theta = sin(2.*angles.Theta());
	
	GDouble phi = angles.Phi();
	
	TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector in lab
	GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
//	Phi = Phi > 0? Phi : Phi + 3.14159;
	
	// polarization BeamProperties
	GDouble Pgamma=polFraction;
	
	if(polAngle == -1)
		Pgamma = 0.;
	
	// SDMEs for 3/2- -> 1/2+ + 0- (doi.org/10.1103/PhysRevC.96.025208)
	GDouble W = 3.*(0.5 - rho011)*sinSqTheta + rho011*(1.+3.*cosSqTheta) - 2.*TMath::Sqrt(3.)*rho031*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*rho03m1*cos(2.*phi)*sinSqTheta;
	
	W -= Pgamma*cos(2.*Phi) * (3.*rho133*sinSqTheta + rho111*(1.+3.*cosSqTheta) - 2.*TMath::Sqrt(3.)*rho131*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*rho13m1*cos(2.*phi)*sinSqTheta);
	
	W -= Pgamma*sin(2.*Phi) * (2.*TMath::Sqrt(3.)*rho231*sin(phi)*sin2Theta + 2.*TMath::Sqrt(3.)*rho23m1*sin(2.*phi)*sinSqTheta);
	
	W *= 1./(4.*PI);
	
// 	return W;
	return complex< GDouble > ( sqrt(fabs(W)) );
}

