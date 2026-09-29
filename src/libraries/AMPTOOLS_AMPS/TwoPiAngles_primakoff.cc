
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <complex.h>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"


#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiAngles_primakoff.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "UTILITIES/BeamProperties.h"

TwoPiAngles_primakoff::TwoPiAngles_primakoff( const vector< string >& args ) :
UserAmplitude< TwoPiAngles_primakoff >( args )
{
	assert( args.size() == 7 );
	
	phipol  = atof(args[0].c_str() )*3.14159/180.; // azimuthal angle of the photon polarization vector in the lab. Convert to radians.
	polFrac  = AmpParameter( args[1] ); // fraction of polarization (0-1)  // use temporarily for beam configuration file

	// BeamProperties configuration file
	// TString beamConfigFile = args[2].c_str();
	// cout << "TwoPiAngles_primakoff  Init: beamConfigFile=" << beamConfigFile << endl;
	TString polFrac_vs_E_fname = args[2].c_str();
	TString polFrac_vs_E_hname = args[3].c_str();
	if (! InitPol && polFrac <= 0) {
          polFraction = 0.; 
	  // cout << "TwoPiAngles_primakoff  Init: polFrac_vs_E_fname=" << polFrac_vs_E_fname << " polFrac_vs_E_hname=" << polFrac_vs_E_hname << endl;
          // TFile* f = new TFile(polFrac_vs_E_fname);
          // polFrac_vs_E = (TH1D*)f->Get(polFrac_vs_E_hname);  
          /*Int_t const nbins=10;
	  Double_t xlo=5.05;
	  Double_t xhi=6.05;
          Double_t xbin_mean[nbins]={5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0};
          Double_t content_mean[nbins]={0.649351,0.672864,0.688222,0.700452,0.717672,0.722537,0.721527,0.658503,0.38531,0.346376}; */
	   
          Int_t const nbins=101;   // Use data from RHO polarization normalized to TPOL
	  Double_t xlo=1.95;
	  Double_t xhi=12.05;
          Double_t xbin_mean[nbins]={2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
                          3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
                          4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
                          5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,
                          6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,
                          7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,
                          8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,
                          9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,
                          10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,
                          11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12};
          Double_t content_mean[nbins]={0,0,0,0,0,0,0,0,0,0.0271335,
                           0.0117895,0.0181829,0.0447745,0.0574396,0.0553967,0.0906021,0.10313,0.144275,0.180627,0.196628,
                           0.244986,0,0,0.186031,0.388328,0.436327,0.483053,0.523823,0.559426,0.593285,
                           0.624764,0.646664,0.667381,0.691165,0.704588,0.71706,0.722791,0.707826,0.545987,0.423112,
                           0.386477,0.298619,0.247285,0.228017,0.226629,0.232654,0.222539,0.238963,0.272686,0.286454,
                           0.329545,0.347663,0.366713,0.39447,0.416947,0.420199,0.438003,0.378953,0.168752,0,
                           0,0,0,0,0,0,0,0,0.0778183,0.0817014,
                           0.101779,0.115131,0.114527,0.0784607,0.0746836,0.0471799,0.0490841,0.0387982,0.0384065,0.0332145,
                           0.0299022,0.0254,0.0141661,0.0243896,0.0231319,0.015943,0.0185516,0.0100371,0.00893803,0,
                           0.0192681,0.00331965,0.0112799,0.0137505,0,0,0,0,0,0,0,}; 
          if (polFrac_vs_E == NULL ) {
	      cout << "TwoPiAngles_primakoff  Init: Create Histogram polFrac_vs_E" << endl;
              polFrac_vs_E = new TH1D ("polFrac_vs_E","TPOL Average Polarizations for CPP",nbins,xlo,xhi);
              polFrac_vs_E->FillN(nbins,xbin_mean,content_mean);
	  }
	  InitPol = true;
	}
	else {
	  cout << "TwoPiAngles_primakoff  Init: Use constant Polarization Fraction=" << polFrac << endl;
	}
	m_rho = atoi( args[4].c_str() );  // Jz component of rho
	PhaseFactor  =atoi( args[5].c_str() ) ;  // prefix factor to amplitudes in computation
	flat = atoi( args[6].c_str() );  // flat=1 uniform angles, flat=0 use YLMs

	cout << "TwoPiAngles_primakoff  Init: phipol=" << phipol << " polFrac=" << polFrac << " m_rho=" << m_rho << " PhaseFactor=" << PhaseFactor << " flat=" << flat << endl;

	assert( ( phipol >= 0.) && (phipol <= 2*3.14159));
	assert( ( polFrac >= 0 ) && ( polFrac <= 1 ) );
        assert( ( m_rho == 1 ) || ( m_rho == 0 ) || ( m_rho == -1 ));
        assert( ( PhaseFactor == 0 ) || ( PhaseFactor == 1 ) || ( PhaseFactor == 2 ) || ( PhaseFactor == 3 ));
	assert( (flat == 0) || (flat == 1) );

	// need to register any free parameters so the framework knows about them
	registerParameter( polFrac );
}


complex< GDouble >
TwoPiAngles_primakoff::calcAmplitude( GDouble** pKin ) const {

	complex< GDouble > i( 0, 1 );
	complex< GDouble > factor( 0, 0 );
	complex< GDouble > Amp( 0, 0 );
	Int_t Mrho=0;

  if (flat == 1) { // no computations needed
     Amp = 1;
     return Amp;
  }


  // for Primakoff, all calculations are in the lab frame. Keep recoil but remember that it cannot be measured by detector.
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector p1     ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p2     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector recoil ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );
	TLorentzVector resonance = p1 + p2;

        TVector3 eps(G_COS(phipol), G_SIN(phipol), 0.0); // beam polarization vector in lab

	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
	TLorentzVector p2_res = resonanceBoost * p2;

        // choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is defined as missing P4
        TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();  
        TVector3 z = -1. * recoil_res.Vect().Unit();
        TVector3 x = y.Cross(z).Unit();
        TVector3 angles1( (p1_res.Vect()).Dot(x),
                         (p1_res.Vect()).Dot(y),
                         (p1_res.Vect()).Dot(z) );
        TVector3 angles2( (p2_res.Vect()).Dot(x),
                         (p2_res.Vect()).Dot(y),
                         (p2_res.Vect()).Dot(z) );


	// Pick pi0 randomly between pi01 and pi02
        TRandom1 *r1 = new TRandom1();
        GDouble CosTheta; 
        if (r1->Rndm() > 0.5) {
	  CosTheta = angles1.CosTheta();
	      }
       else {
	  CosTheta = angles2.CosTheta();
	      }
   
        GDouble phi = angles1.Phi();
        // GDouble sinSqTheta = G_SIN(angles.Theta())*G_SIN(angles.Theta());
        // GDouble sin2Theta = G_SIN(2.*angles.Theta());

        GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

        GDouble psi = Phi - phi;               // define angle difference 
        if(psi < -1*PI) psi += 2*PI;
        if  (psi > PI) psi -= 2*PI;

	/*cout << " recoil_res Angles="; recoil_res.Vect().Print();
	cout << " p1_res Angles="; p1_res.Vect().Print();
	cout << "phi= " << phi << endl;
	cout << " psi=" << psi << endl;*/
	double polFraction = 0.;  // Get energy-dependenet polarization fraction

	if (polFrac <= 0) {
           int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
           if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
	     polFraction = 0.;
           }
	   else {
	      polFraction = polFrac_vs_E->GetBinContent(bin);
            }
	  }
	  else {
	    polFraction = polFrac;
	  }
	// cout << " beam E=" << beam.E() <<  " polFrac =" << polFrac << " polFraction=" << polFraction << endl;

	switch (PhaseFactor) {
        case 0:
	  Mrho = m_rho;
	  Amp = G_SQRT(1-polFraction)*(-G_SIN(Phi)* Y( 0, Mrho, CosTheta, phi) );
	  break;
        case 1:
	  Mrho = m_rho;
	  Amp = G_SQRT(1+polFraction)*(G_COS(Phi)* Y( 0, Mrho, CosTheta, phi)  );
	  break;
        case 2:
	  Mrho = m_rho;
	  factor = exp(-i*Phi)* Y( 1, Mrho, CosTheta, phi);
	  Amp = G_SQRT(1-polFraction)* imag(factor);
	  break;
        case 3:
	  Mrho = m_rho;
	  factor = exp(-i*Phi)* Y( 1, Mrho, CosTheta, phi);
	  Amp = G_SQRT(1+polFraction)* real(factor);
	  break;
	}


	if (abs(Amp) <= 0) cout << " m_rho=" << m_rho << " CosTheta=" << CosTheta << " phi=" << phi << " factor=" << factor << " Amp=" << Amp << endl;

	return Amp;
}

