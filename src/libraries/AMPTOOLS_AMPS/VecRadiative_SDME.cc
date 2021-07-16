#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
//#include "AMPTOOLS_AMPS/VecRadiative_SDME.h"
//#include "UTILITIES/BeamProperties.h"

#include "GlueXAmp/VecRadiative_SDME.h"

VecRadiative_SDME::VecRadiative_SDME( const vector< string >& args ) :
UserAmplitude< VecRadiative_SDME >( args )
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
}

complex< GDouble >
VecRadiative_SDME::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {
  
  GDouble sinSqTheta = userVars[kSinSqTheta];
  GDouble sin2Theta = userVars[kSin2Theta];
  GDouble cosTheta = userVars[kCosTheta];
  GDouble phi = userVars[kPhi];
  GDouble bigPhi = userVars[kBigPhi];
  GDouble polFrac = userVars[kPolFrac];
  
  GDouble m_rho110 = 0.5 * ( 1 - m_rho000 );
  
  GDouble W = 1 - sinSqTheta * m_rho110 - cosTheta*cosTheta*m_rho000 +
              sinSqTheta*cos(2*phi)*m_rho1m10 +
              sqrt(2)*m_rho100*sin2Theta*cos(phi);
  
  W -= polFrac * cos(2*bigPhi) * ( 2*m_rho111 + sinSqTheta*(m_rho001-m_rho111) +
                                   sinSqTheta*cos(2*phi)*m_rho1m11 +
                                   sqrt(2)*m_rho101*sin2Theta*cos(phi) );
  
  W += polFrac * sin(2*bigPhi) * ( m_rho1m12*sinSqTheta*sin(2*phi) +
                                   sqrt(2)*m_rho102*sin2Theta*sin(phi) );
  
  W *= 3/(8*PI);
  
  return complex< GDouble >( sqrt(fabs(W)) );
}


void
VecRadiative_SDME::calcUserVars( GDouble** pKin, GDouble* userVars ) const {
  
  // written in the context of omega photoproduction with omega -> gamma  pi0,
  // but these angle definitions in the helicity frame for a two-body decay
  // are fairly standard
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
  TLorentzVector pi0    ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
  TLorentzVector gamma  ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );
  
  TLorentzVector omega = pi0 + gamma;
  
  // construct a boost that will take us to the omega rest frame:
  TLorentzRotation omegaRestBoost( -omega.BoostVector() );
  
  TLorentzVector recoil_omegaRest = omegaRestBoost * recoil;
  TLorentzVector gamma_omegaRest = omegaRestBoost * gamma;
  
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
  
  // a unit vector for the polarization direction in the lab
  TVector3 polUnitVec( cos( m_polAngle*TMath::DegToRad() ),
                       sin( m_polAngle*TMath::DegToRad() ), 0.0 );
  
  // now get the angle between the production plane (defined by the
  // y axis in the helicity frame) and the photon polarization vector:
 
  GDouble bigPhi = atan2( yHel.Dot( polUnitVec ),
                          beam.Vect().Unit().Dot( polUnitVec.Cross( yHel ) ) );

  // get the angles of the radiated photon in the helicity frame:
  GDouble cosTheta = gamma_omegaRest.Vect().Unit().Dot( zHel );
  GDouble theta = acos( cosTheta );
  
  GDouble cosPhi =  yHel.Dot( zHel.Cross( gamma_omegaRest.Vect() ).Unit() );
  GDouble sinPhi = -xHel.Dot( zHel.Cross( gamma_omegaRest.Vect() ).Unit() );
  GDouble phi = atan2( sinPhi, cosPhi );

  userVars[kCosTheta]   = cosTheta;
  userVars[kSinSqTheta] = sin(theta)*sin(theta);
  userVars[kSin2Theta]  = sin(2*theta);
  userVars[kPhi] = phi;
  userVars[kBigPhi] = bigPhi;
  
  GDouble polFrac;
  if(m_polFraction > 0.) { // for fitting with constant polarization

    polFrac = m_polFraction;
  }
  else{

    int bin = m_polFrac_vs_E->GetXaxis()->FindBin( pKin[0][0] );
    if( bin == 0 || bin > m_polFrac_vs_E->GetXaxis()->GetNbins() ) polFrac = 0.;
    else polFrac = m_polFrac_vs_E->GetBinContent( bin );
  }
  
  userVars[kPolFrac] = polFrac;
}

#ifdef GPU_ACCELERATION
void
VecRadiative_SDME::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  GPUVecRadiative_SDME_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
                             m_rho000, m_rho100, m_rho1m10,
                             m_rho111, m_rho001, m_rho101,
                             m_rho1m11, m_rho102, m_rho1m12 );
}
#endif

