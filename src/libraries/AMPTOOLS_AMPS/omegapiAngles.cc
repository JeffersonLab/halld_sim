#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <cmath>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TMath.h"
#include "omegapiAngles.h"


vector <double> getVectorDecayAngles( const TLorentzVector& beamPLab,
                                      const TLorentzVector& particleXLab,
                                      const TLorentzVector& parentLab, 
                                      const TLorentzVector& firstDaughterLab, 
                                      const TLorentzVector& secondDaughterLab){
  // Calculating the angles of the normal to the piplus+piminus plane 
  // beamPLab = beam + proton (center of mass), particleXLab = b1, 
  // parentLab = omega, firstDaughterLab = piplus, secondDaughterLab = piminus

  // Boost all particles to the (beamPLab) center-of-mass frame (cm)
  TLorentzVector particleX_cm = particleXLab;
  TLorentzVector parent_cm = parentLab;
  TLorentzVector firstDaughter_cm = firstDaughterLab;
  TLorentzVector secondDaughter_cm = secondDaughterLab;
  TVector3 cmBoostVector = -beamPLab.BoostVector();
  particleX_cm.Boost(cmBoostVector);
  parent_cm.Boost(cmBoostVector);
  firstDaughter_cm.Boost(cmBoostVector);
  secondDaughter_cm.Boost(cmBoostVector);

  // Boost parent and daughters to rest frame of particle X (x)
  TLorentzVector parent_x = parent_cm;
  TLorentzVector firstDaughter_x = firstDaughter_cm;
  TLorentzVector secondDaughter_x = secondDaughter_cm;
  TVector3 xBoost = -particleX_cm.BoostVector();
  parent_x.Boost(xBoost);
  firstDaughter_x.Boost(xBoost);
  secondDaughter_x.Boost(xBoost);
  
  // Boost daughters to parent's rest frame (p)
  TLorentzVector firstDaughter_p = firstDaughter_x;
  TLorentzVector secondDaughter_p = secondDaughter_x;
  TVector3 parentBoost = -parent_x.BoostVector();
  firstDaughter_p.Boost(parentBoost);
  secondDaughter_p.Boost(parentBoost);

  // Define the unit vectors for the helicity frame of the parent
  TVector3 firstDaughter_p_unit = (firstDaughter_p.Vect()).Unit();
  TVector3 secondDaughter_p_unit = (secondDaughter_p.Vect()).Unit();
  TVector3 normal_p_unit = firstDaughter_p_unit.Cross(secondDaughter_p_unit);


  // The axis of this helicity frame are defined as follows:
  // z-axis is along the direction of the parent in the X rest frame
  // y-axis is normal to the decay plane defined by the parent and particle X
  //        the parent's vector is in the X rest frame and the particle X's
  //        vector is in the cm frame 
  TVector3 z = (parent_x.Vect()).Unit();
  TVector3 y = ((particleX_cm.Vect()).Cross(z)).Unit();
  TVector3 x = (y.Cross(z)).Unit();

  // The decay vector is normal to decay plane for omega->3pi
  // In the case of vec->ps1+ps2, the decay vector is one of the ps
  TVector3 decayVector;
  if(secondDaughterLab.E() > 0) decayVector = normal_p_unit;
  else decayVector = firstDaughter_p_unit;
  
  TVector3 components(decayVector.Dot(x),decayVector.Dot(y),decayVector.Dot(z));

  double thetaHelicity = components.Theta();
  double phiHelicity = components.Phi();

  // Compute the variable lambda for the omega->3pi decay
  TVector3 daughterCross = (firstDaughter_p.Vect()).Cross(secondDaughter_p.Vect());
  double m0 = 0.1349766; // mass of pi0
  double mq = 0.1395702; // mass of charged pions
  double threePiMass = parentLab.M2();
  double lambda_max =  (3.0/4.0) *
      TMath::Power( (1.0/9.0) * (5*threePiMass + 3*(m0*m0 - 4*mq*mq) -
                   4*sqrt(threePiMass * threePiMass +
                   3*threePiMass * (m0*m0 - mq*mq))), 2);
  double lambda = fabs(daughterCross.Mag2()) / lambda_max;

  return {thetaHelicity, phiHelicity, lambda};

}
vector<double> getXDecayAngles( double polAngle, 
                                const TLorentzVector& beamLab, 
                                const TLorentzVector& beamPLab,
                                const TLorentzVector& particleXLab,
                                const TLorentzVector& daughterLab){
    // Calculating the helicity angles of omega in  the b1 decay
    // daughterLab = omega, particleXLab = b1, beamLab = beam,
    // beamPLab = beam + proton (center of mass)

    // Boost all particles to the (beamPLab) center-of-mass frame (cm)
    TLorentzVector beam_cm = beamLab;
    TLorentzVector particleX_cm = particleXLab;
    TLorentzVector daughter_cm = daughterLab;
    TVector3 cmBoostVector = -beamPLab.BoostVector();
    beam_cm.Boost(cmBoostVector);
    particleX_cm.Boost(cmBoostVector);
    daughter_cm.Boost(cmBoostVector);

    // Boost daughter to particle X rest frame (x)
    TLorentzVector daughter_x = daughter_cm;
    TVector3 xBoost = -particleX_cm.BoostVector();
    daughter_x.Boost(xBoost);

    // get the unit vectors in space
    TVector3 daughter_x_unit = (daughter_x.Vect()).Unit();
    TVector3 beam_cm_unit = (beam_cm.Vect()).Unit();

    // The axis of this helicity frame are defined as follows:
    // z-axis is along the direction of particle X in the center-of-mass frame
    // y-axis is normal to the production plane defined by particle X and
    //        the beam both vectors are expressed in the center-of-mass frame
    TVector3 z = (particleX_cm.Vect()).Unit();
    TVector3 y = ((beam_cm_unit).Cross(z)).Unit();
    TVector3 x = (y.Cross(z)).Unit();

    // One could use the Gottfried-Jackson frame instead of the helicity frame
    // The axis would be defined as follows:
    // z-axis is along the direction of the beam in the particle X rest frame
    // y-axis is normal to the production plane defined by particle X and
    //        the beam both vectors are expressed in the center-of-mass frame
    //        (the same as helicity frame)
    // TLorentzVector beam_x = beam_cm;
    // beam_x.Boost(-1.0*xBoost);
    // TVector3 z = (beam_x.Vect()).Unit();
    // TVector3 y = ((beam_cm_unit.Vect()).Cross(particleX_cm.Vect())).Unit();
    // TVector3 x = (y.Cross(z)).Unit();

    TVector3 components(daughter_x_unit.Dot(x), daughter_x_unit.Dot(y),
                        daughter_x_unit.Dot(z));

    double theta = components.Theta();
    double phi = components.Phi();

    // Compute the production angle (bigPhi) between the polarization 
    // angle and the normal to the production plane
    // But first, make sure the polarization angle is in radians
    static bool warned = false; // only warn once
    if(!warned && (fabs(polAngle) > 2 * TMath::Pi())){
      cerr << "[Notice] getXDecayAngles(): polAngle = " << polAngle
           << " appears to be in degrees. Converting to radians."
           << endl;
      polAngle = DEG_TO_RAD * polAngle;
      warned = true;
    }
    TVector3 eps(cos(polAngle), sin(polAngle), 0.0); // polarization vector
    double bigPhi = atan2(y.Dot(eps), beamLab.Vect().Unit().Dot(eps.Cross(y)));    

    return {theta, phi, bigPhi};
}
