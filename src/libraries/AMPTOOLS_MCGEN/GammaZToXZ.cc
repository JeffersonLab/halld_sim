/*
 *  GammaZToXZ.cc
 *  GlueXTools
 *
 * Copied from GammaPToXP.    Elton 5/26/2020
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 * Copied from 
 */

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

#include "particleType.h"

#include "AMPTOOLS_MCGEN/GammaZToXZ.h"
#include "UTILITIES/BeamProperties.h"

GammaZToXZ::GammaZToXZ( float massX, TString beamConfigFile, Double_t Bslope) :
m_target( 0, 0, 0, ParticleMass(Pb208)),
m_slope( Bslope ),
m_recoil (ParticleMass(Pb208)),
m_childMass( 0 )
{
  kMproton=ParticleMass(Proton);
  kMneutron=ParticleMass(Neutron);
  // kMZ = 108.;      //  mass of Sn116 
  // kMZ = 208.*0.931494;      //  use mass of Pb as it is in the particle table
  kMZ = ParticleMass(Pb208);      //  use mass of Pb as it is in the particle table
  kMPion = ParticleMass(PiPlus);
  kMPi0 = ParticleMass(Pi0);
  kMKaon = ParticleMass(KPlus);

  m_childMass.push_back( massX );

  // get beam properties from configuration file
  BeamProperties beamProp(beamConfigFile);
  cobrem_vs_E = (TH1D*)beamProp.GetFlux();
  cobrem_vs_E->GetName();
}

Kinematics* 
GammaZToXZ::generate(){


  const double kPi = 3.14159;

// Double_t kMProton = ParticleMass(Proton);
  double EtaMass=m_childMass[0];
  double recoilMass = m_recoil;
// Double_t masses[2] = {EtaMass,recoilMass};

  double beamE = cobrem_vs_E->GetRandom();
  TLorentzVector beam;
  beam.SetPxPyPzE(0,0,beamE,beamE);
  TLorentzVector cm = beam + m_target;
  
  TLorentzRotation lab2cmBoost( -cm.BoostVector() );
  TLorentzRotation cm2labBoost( cm.BoostVector() );

   TLorentzVector cmCM = lab2cmBoost * cm;
  
  double cmEnergy = ( lab2cmBoost * cm ).E();
  double beamMomCM = cmMomentum( cmEnergy, beam.M(), m_target.M() );
  double EtaMomCM  = cmMomentum( cmEnergy, EtaMass, recoilMass );

/*cout  << "M=" << beam.M() << " beam="; beam.Print();
  cout << "M=" << m_target.M() << " m_target="; m_target.Print();
  cout << "M=" << cm.M() << " cm="; cm.Print();*/

  // generate an exponetial t-slope between t1 and t2

  double t, tMaxkin, tMin, tMax;
    // generate the t-distribution. t is positive here (i.e. should be -t)
  
    tMaxkin = 4. * beamMomCM * EtaMomCM;
    tMax = 0.2;   // restrict max to make more efficient for Primakoff generation (about 2. deg at 0.05 GeV-2)
    // tMax = 1.;   // restrict max to make more efficient for Primakoff generation
    tMin = abs(pow(EtaMass,4)/(2*cmEnergy*cmEnergy) - (beamMomCM-EtaMomCM)*(beamMomCM-EtaMomCM));  // treating t as positive

    double Irandom = gRandom->Uniform();
    // generate random t with exponential between tMin and tMax
    t = -log( Irandom*(exp(-m_slope*tMax) - exp(-m_slope*tMin)) + exp(-m_slope*tMin))/m_slope;

   TVector3 EtaMom3CM;
   double thetaCM = 2.*sqrt((t- tMin)/tMaxkin); // acos( 1. - 2.*t/tMax ) -> use small angle approximation to avoid roundoff. For heavy target, CM=Lab. Opposite to target
	// double thetaCM = acos( 1. - 2.*t/tMaxkin );
   double phiCM = random( -kPi, kPi ); 

   EtaMom3CM.SetMagThetaPhi( EtaMomCM, thetaCM, phiCM);
	
   TLorentzVector EtaMom4CM( EtaMom3CM, sqrt( EtaMom3CM.Mag2() + EtaMass * EtaMass ) );
   TLorentzVector recoilMom4CM = cmCM - EtaMom4CM;


   
   TLorentzVector EtaMom4Lab = cm2labBoost * EtaMom4CM;
   TLorentzVector recoil = cm2labBoost * recoilMom4CM;

double trecoil = 2*recoil.M()*(recoil.E()-recoil.M());

// cout << " tMin=" << tMin << " tMax=" << tMax << " t=" << t << " trecoil=" << trecoil << " m_slope=" << m_slope << " thetaCM=" << thetaCM*180/kPi << " phiCM=" << phiCM*180./kPi << endl;

  vector< TLorentzVector > allPart;
  allPart.push_back( beam );
  allPart.push_back( EtaMom4Lab );
  allPart.push_back( recoil );

/*cout  << "M=" << beam.M() << " beam="; beam.Print();
  cout  << "cmCM=" << cmCM.M() << " cmCM="; cmCM.Print();
  cout << "M=" << recoil.M() << " recoil="; recoil.Print();
  cout << " m_slope=" << m_slope << " M=" << EtaMom4Lab.M() << " EtaMom4Lab="; EtaMom4Lab.Print(); cout << endl;*/

  double genWeight = 1.;
 
  return new Kinematics( allPart, genWeight );
}
double
GammaZToXZ::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
GammaZToXZ::random( double low, double hi ) const {

        return( ( hi - low ) * gRandom->Uniform() + low );
}

