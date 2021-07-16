
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
//#include "AMPTOOLS_AMPS/ThreePiAnglesSchilling.h"
//#include "AMPTOOLS_AMPS/clebschGordan.h"
//#include "AMPTOOLS_AMPS/wignerD.h"

#include "ThreePiAnglesSchilling.h"

//#include "UTILITIES/BeamProperties.h"

ThreePiAnglesSchilling::ThreePiAnglesSchilling( const vector< string >& args ) :
UserAmplitude< ThreePiAnglesSchilling >( args )
{
  assert( args.size() == 11 );
  
  m_rho000  = AmpParameter( args[0] );
  m_rho100  = AmpParameter( args[1] );
  m_rho1m10 = AmpParameter( args[2] );
  
  m_rho111  = AmpParameter( args[3] );
  m_rho001  = AmpParameter( args[4] );
  m_rho101  = AmpParameter( args[5] );
  m_rho1m11 = AmpParameter( args[6] );
  
  m_rho102  = AmpParameter( args[7] );
  m_rho1m12 = AmpParameter( args[8] );
  
  m_polAngle = AmpParameter( args[9] );
  
  m_polFraction = atof(args[10].c_str());
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_rho000 );
  registerParameter( m_rho100 );
  registerParameter( m_rho1m10 );
  
  registerParameter( m_rho111 );
  registerParameter( m_rho001 );
  registerParameter( m_rho101 );
  registerParameter( m_rho1m11 );
  
  registerParameter( m_rho102 );
  registerParameter( m_rho1m12 );
  
  registerParameter( m_polAngle );
  /*
   if (m_polFraction > 0.0)
   cout << "Fitting with constant polarization" << endl;
   else
   {
   cout << "Fitting with polarization from BeamProperties class" << endl;
   // BeamProperties configuration file
   TString beamConfigFile = args[10].c_str();
   BeamProperties beamProp(beamConfigFile);
   m_polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
   }
   */
  
  TFile* f = new TFile( "TPol_201808.root" );
  m_polFrac_vs_E = (TH1D*)f->Get( "hPol0" );
}

complex< GDouble >
ThreePiAnglesSchilling::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  GDouble sinSqTheta = userVars[kSinSqTheta];
  GDouble sin2Theta = userVars[kSin2Theta];
  GDouble cosTheta = userVars[kCosTheta];
  GDouble phi = userVars[kPhi];
  GDouble bigPhi = userVars[kBigPhi];
  GDouble polFrac = userVars[kPolFrac];
  
  // vector meson production from K. Schilling et. al.

  GDouble W = 0.5*(1 - m_rho000) + 0.5*(3*m_rho000 - 1)*cosTheta*cosTheta -
              sqrt(2)*m_rho100*sin2Theta*cos(phi) - m_rho1m10*sinSqTheta*cos(2*phi);

  W -= polFrac * cos(2*bigPhi) * ( m_rho111*sinSqTheta + m_rho001*cosTheta*cosTheta -
                                   sqrt(2)*m_rho101*sin2Theta*cos(phi) -
                                   m_rho1m11*sinSqTheta*cos(2*phi) );

  W -= polFrac * sin(2*bigPhi) * ( sqrt(2)*m_rho102*sin2Theta*sin(phi) +
                                   m_rho1m12*sinSqTheta*sin(2*phi));
  W *= 3/(4*PI);

  return complex< GDouble > ( sqrt(fabs(W)) );
}

void
ThreePiAnglesSchilling::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );
  TLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] );
  
  TLorentzVector omega = p1 + p2 + p3;
  
  // construct a boost that will take us to the omega rest frame:
  TLorentzRotation omegaRestBoost( -omega.BoostVector() );
  
  TLorentzVector recoil_omegaRest = omegaRestBoost * recoil;
  TLorentzVector p1_omegaRest = omegaRestBoost * p1;
  TLorentzVector p2_omegaRest = omegaRestBoost * p2;

  // for the 3pi decay the normal to the decay plane is the
  // key vector that describes the decay
  TVector3 norm = (p2_omegaRest.Vect().Cross(p1_omegaRest.Vect())).Unit();
  
  // in the helicity frame the z axis is opposite the recoil in the CM
  // frame, but the direction of the recoil in the CM frame is the same
  // as the direction of the recoil in the omega rest frame, so this
  // works fine:
  TVector3 zHel = -1*recoil_omegaRest.Vect().Unit();
 
  // the y axis is normal to the production plane production plane, which
  // can be obtained in a variety of ways and is unchanged by boosts
  // to the CM or omega rest frame as all of these boost vectors
  // lie in the production plane
  TVector3 yHel = (beam.Vect().Cross(omega.Vect())).Unit();
 
  // and the x axis is constructed so the coordinate system is
  // right handed:
  TVector3 xHel = yHel.Cross(zHel);
 
  // construct the normal to the decay plane in the helicity frame:
  TVector3 normHel( norm.Dot( xHel ), norm.Dot( yHel ), norm.Dot( zHel ) );
  
  userVars[kCosTheta] = normHel.CosTheta();
  userVars[kSinSqTheta] = sin(normHel.Theta())*sin(normHel.Theta());
  userVars[kSin2Theta] = sin(2*normHel.Theta());
  userVars[kPhi] = normHel.Phi();
  
  // a unit vector for the polarization direction in the lab
  TVector3 polUnitVec( cos( m_polAngle*TMath::DegToRad() ),
                       sin( m_polAngle*TMath::DegToRad() ), 0.0 );
  
  // now get the angle between the production plane (defined by the
  // y axis in the helicity frame) and the photon polarization vector:
  userVars[kBigPhi] = atan2( yHel.Dot( polUnitVec ),
                             beam.Vect().Unit().Dot( polUnitVec.Cross( yHel ) ) );

  GDouble polFrac = 0;
  if( m_polFraction > 0 ) {

    // for fitting with constant polarization
    polFrac = m_polFraction;
  }
  else{
    int bin = m_polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
    if( bin == 0 || bin > m_polFrac_vs_E->GetXaxis()->GetNbins() ) polFrac = 0;
    else polFrac = m_polFrac_vs_E->GetBinContent(bin);
  }
  
  userVars[kPolFrac] = polFrac;
 }

