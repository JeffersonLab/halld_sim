#if !defined(GAMMAZTOXZ)
#define GAMMAZTOXZ

/*
 *  GammaZToXP.h
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 *  Copied from GAmmaPToXP.h    5/26/2020
 */

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"

#include "IUAmpTools/Kinematics.h"

class Kinematics;

class GammaZToXZ {
  
public:

  enum Recoil { kProton, kNeutron, kZ};
  
  GammaZToXZ( float massX, TString beamConfigFile, Double_t Bslope);
  
  Kinematics* generate();
  
private:
  
  double kMproton,kMneutron,kMZ, kMPion, kMKaon, kMPi0;
  TLorentzVector m_beam;
  TLorentzVector m_target;
  double m_slope;
  double m_recoil;
  
  vector< double > m_childMass;
  
  TH1D *cobrem_vs_E;
  
  double cmMomentum( double M, double m1, double m2 ) const;
  double random( double low, double hi ) const;

};

#endif
