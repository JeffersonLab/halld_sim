#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/DeltaAngles_Rho.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

DeltaAngles_Rho::DeltaAngles_Rho( const vector< string >& args ) :
UserAmplitude< DeltaAngles_Rho >( args )
{
	assert( args.size() == 12 || args.size() == 11 );
	
	delta_rho011  = AmpParameter( args[0] );
	delta_rho031  = AmpParameter( args[1] );
	delta_rho03m1 = AmpParameter( args[2] );
	
	delta_rho111  = AmpParameter( args[3] );
	delta_rho133  = AmpParameter( args[4] );
	delta_rho131  = AmpParameter( args[5] );
	delta_rho13m1 = AmpParameter( args[6] );
	
	delta_rho231  = AmpParameter( args[7] );
	delta_rho23m1 = AmpParameter( args[8] );

	frame = string( args[9] ); // choose between Gottfried Jackson (GJ) and helicity (Hel) frame
	
	// need to register any free parameters so the framework knows about them
	registerParameter( delta_rho011 );
	registerParameter( delta_rho031 );
	registerParameter( delta_rho03m1 );
	
	registerParameter( delta_rho111 );
	registerParameter( delta_rho133 );
	registerParameter( delta_rho131 );
	registerParameter( delta_rho13m1 );
	
	registerParameter( delta_rho231 );
	registerParameter( delta_rho23m1 );

	std::cout<< "number of arguments:" <<args.size()<< std::endl;
	
	if(args.size() == 12){
		polAngle  = atof(args[10].c_str() ); // azimuthal angle of the photon polarization vector in the lab.
		polFraction = AmpParameter( args[11] ); // fraction of polarization (0-1)
		std::cout << "Fixed polarisation of " << polFraction << " and angle of " << polAngle << " degrees." << std::endl;
	}
	else
		assert(0);
}

complex< GDouble >
DeltaAngles_Rho::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {	

	GDouble sinSqTheta 	= userVars[kSinSqTheta];
//	GDouble cosSqTheta 	= userVars[kCosSqTheta];
	GDouble sin2Theta	= userVars[kSin2Theta];
	GDouble phi			= userVars[kPhi];
	GDouble cosTheta	= userVars[kCosTheta];
	GDouble bigPhi		= userVars[kBigPhi];
	GDouble Pgamma		= userVars[kPgamma];

	//SDMEs for 3/2- -> 1/2+ + 0- (doi.org/10.1103/PhysRevC.96.025208)

	GDouble W = 3.*(0.5 - delta_rho011)*sinSqTheta + delta_rho011*(1.+3.*cosTheta*cosTheta) - 2.*TMath::Sqrt(3.)*delta_rho031*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*delta_rho03m1*cos(2.*phi)*sinSqTheta;

	W -= Pgamma*cos(2.*bigPhi) * (3.*delta_rho133*sinSqTheta + delta_rho111*(1.+3.*cosTheta*cosTheta) - 2.*TMath::Sqrt(3.)*delta_rho131*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*delta_rho13m1*cos(2.*phi)*sinSqTheta);

	W -= Pgamma*sin(2.*bigPhi) * (2.*TMath::Sqrt(3.)*delta_rho231*sin(phi)*sin2Theta + 2.*TMath::Sqrt(3.)*delta_rho23m1*sin(2.*phi)*sinSqTheta);

	W *= 1./(4.*PI);

	

// 	return W;
	return complex< GDouble > ( sqrt(fabs(W)) );
}

void
DeltaAngles_Rho::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
	
	TLorentzVector target ( 0, 0, 0, 0.9382720813);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  	TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] ); 
	
  	TLorentzVector recoil = p2 + p3;

	TLorentzVector resonance = proton + p1;
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector target_res = resonanceBoost * target;
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
	
	// normal to the production plane
	// choose Gottfried-Jackson frame: z-axis along target direction in baryon rest frame
	TVector3 y = (beam_res.Vect().Unit().Cross(recoil_res.Vect().Unit())).Unit();
	TVector3 z = target_res.Vect().Unit();

	// TVector3 y = (target_res.Vect().Unit().Cross(recoil_res.Vect().Unit())).Unit();
	// TVector3 z = target_res.Vect().Unit();
	TVector3 x = y.Cross(z).Unit();
	
	TVector3 angles( (p1_res.Vect()).Dot(x),
					(p1_res.Vect()).Dot(y),
					(p1_res.Vect()).Dot(z) );

	GDouble phi = angles.Phi();
	GDouble cosTheta = angles.CosTheta();
	
	GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
	// GDouble cosSqTheta = cos(angles.Theta())*cos(angles.Theta());
	GDouble sin2Theta = sin(2.*angles.Theta());

	
	TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector in lab
	GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	//Phi = Phi > 0? Phi : Phi + 3.14159;
	
	// polarization BeamProperties
	GDouble Pgamma=polFraction;
	
	if(polAngle == -1)
		Pgamma = 0.;

	//angles in the helicity frame:
	TVector3 z_hel = -recoil_res.Vect().Unit();
	TVector3 y_hel = (beam_res.Vect().Cross(recoil_res.Vect())).Unit();
	TVector3 x_hel = (y_hel.Cross(z_hel)).Unit();

	TVector3 p1_hel(p1_res.Vect().Dot(x_hel),p1_res.Vect().Dot(y_hel),p1_res.Vect().Dot(z_hel));

	GDouble phi_hel = p1_hel.Phi();
	GDouble cosTheta_hel = p1_hel.CosTheta();

	GDouble sinSqTheta_hel = sin(p1_hel.Theta())*sin(p1_hel.Theta());
	// GDouble cosSqTheta_hel = cos(p2_hel.Theta())*cos(p2_hel.Theta());
	GDouble sin2Theta_hel = sin(2.*p1_hel.Theta());

	if(frame=="Hel"){
		userVars[kCosTheta]	= cosTheta_hel;
		userVars[kSinSqTheta] = sinSqTheta_hel;
		userVars[kSin2Theta] = sin2Theta_hel;
		userVars[kPhi] = phi_hel;
	}
	else{
		userVars[kCosTheta]	= cosTheta;
		userVars[kSinSqTheta] = sinSqTheta;
		userVars[kSin2Theta] = sin2Theta;
		userVars[kPhi] = phi;
	} 

	userVars[kPgamma] = Pgamma;
	userVars[kBigPhi] = Phi;
	
}

#ifdef GPU_ACCELERATION
void
DeltaAngles_Rho::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

	GPUDeltaAngles_Rho_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
		delta_rho011, delta_rho031, delta_rho03m1, delta_rho111, delta_rho133, delta_rho131, delta_rho13m1, delta_rho231, delta_rho23m1, polAngle );
}

#endif // GPU_ACCELERATION
