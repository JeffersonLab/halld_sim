
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Lambda1520Angles.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

Lambda1520Angles::Lambda1520Angles( const vector< string >& args ) :
UserAmplitude< Lambda1520Angles >( args )
{
	assert( args.size() == 9 ||  args.size() == 11 || args.size() == 13 );
	
	rho011  = AmpParameter( args[0] );
	rho031  = AmpParameter( args[1] );
	rho03m1 = AmpParameter( args[2] );
	
	rho111  = AmpParameter( args[3] );
	rho133  = AmpParameter( args[4] );
	rho131  = AmpParameter( args[5] );
	rho13m1 = AmpParameter( args[6] );
	
	rho231  = AmpParameter( args[7] );
	rho23m1 = AmpParameter( args[8] );
	
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
	
	if( args.size() == 9 ) {
		// 1. polarization information must be included in beam photon four vector
		//    Usage: amplitude <reaction>::<sum>::<ampName>

		polInTree = true;
	} else if( args.size() == 11 ) {
		// 2. polarization fixed per amplitude and passed as flag
		//    Usage: amplitude <reaction>::<sum>::<ampName> <polAngle> <polFraction>
		polInTree = false;
		polAngle = atof( args[9].c_str() );
		polFraction = atof( args[10].c_str() );
	} else {
		// 2. polarization fixed per amplitude and passed as flag
		//    Usage: amplitude <reaction>::<sum>::<ampName> <polAngle> <polFraction=0.> <rootFile> <hist>
		polInTree = false;
		polAngle = atof( args[9].c_str() );	
		polFraction = 0.;
		TFile* f = new TFile( args[11].c_str() );
        polFrac_vs_E = (TH1D*)f->Get( args[12].c_str() );
        assert( polFrac_vs_E != NULL );
	}
}


complex< GDouble >
Lambda1520Angles::calcAmplitude( GDouble** pKin ) const {
	
	TLorentzVector target ( 0, 0, 0, 0.9382720813);

	TLorentzVector beam;
    TVector3 eps;
   	if(polInTree) {
    	beam.SetPxPyPzE( 0., 0., pKin[0][0], pKin[0][0]);
    	eps.SetXYZ(pKin[0][1], pKin[0][2], 0.); // makes default output gen_amp trees readable as well (without transforming)
   	} else {
    	beam.SetPxPyPzE( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
    	eps.SetXYZ(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
   	}

	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	
	TLorentzVector resonance = p1 + p2;
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector target_res = resonanceBoost * target;
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
	
	// normal to the production plane
	TVector3 y = (beam_res.Vect().Unit().Cross(recoil_res.Vect().Unit())).Unit();
	// choose Gottfried-Jackson frame: z-axis along -target direction in baryon rest frame
	TVector3 z = -1. * target_res.Vect().Unit();
	TVector3 x = y.Cross(z).Unit();
	
	TVector3 angles( (p1_res.Vect()).Dot(x),
					(p1_res.Vect()).Dot(y),
					(p1_res.Vect()).Dot(z) );
	
	GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
	GDouble cosSqTheta = cos(angles.Theta())*cos(angles.Theta());
	GDouble sin2Theta = sin(2.*angles.Theta());
	
	GDouble phi = angles.Phi();
	
	GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));
	Phi = Phi > 0? Phi : Phi + 3.14159;
	
	// get beam polarization
	GDouble Pgamma;
	if(polInTree) {
		Pgamma = eps.Mag();
	} else {
		if(polFraction > 0.) { // for fitting with constant polarization 
			Pgamma = polFraction;
		} else{
			int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
			if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
				Pgamma = 0.;
			} else 
				Pgamma = polFrac_vs_E->GetBinContent(bin);
		}
	}

	
	// SDMEs for 3/2- -> 1/2+ + 0- (doi.org/10.1103/PhysRevC.96.025208)
	GDouble W = 3.*(0.5 - rho011)*sinSqTheta + rho011*(1.+3.*cosSqTheta) - 2.*TMath::Sqrt(3.)*rho031*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*rho03m1*cos(2.*phi)*sinSqTheta;
	
	W -= Pgamma*cos(2.*Phi) * (3.*rho133*sinSqTheta + rho111*(1.+3.*cosSqTheta) - 2.*TMath::Sqrt(3.)*rho131*cos(phi)*sin2Theta - 2.*TMath::Sqrt(3.)*rho13m1*cos(2.*phi)*sinSqTheta);
	
	W -= Pgamma*sin(2.*Phi) * (2.*TMath::Sqrt(3.)*rho231*sin(phi)*sin2Theta + 2.*TMath::Sqrt(3.)*rho23m1*sin(2.*phi)*sinSqTheta);
	
	W *= 1./(4.*PI);
	
// 	return W;
	return complex< GDouble > ( sqrt(fabs(W)) );
}

