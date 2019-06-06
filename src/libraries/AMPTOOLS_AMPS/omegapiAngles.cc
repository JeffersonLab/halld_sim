//The purpose of this function is to return the decay angles of the daughter particle in the helicity frame of the parent.
//For a two particle decay, the function arguments are the lab frame TLorentzVectors of: the daughter particle, the parent particle, the inverse of the X-axis and the reference Z-axis.
//For a three particle decay, the function arguments are the same with the addition of the second daughter particle.
//Please note that in the three particle decay case the normal to the decay plane is calculated as the cross product of momentum vector of the first daughter and the second daughter
//In theory the function could also return the decay angles in other frames when called with the appropriate argument change.
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

#include "omegapiAngles.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

vector <double> getomegapiAngles(TLorentzVector daughter, TLorentzVector parent, TLorentzVector InverseOfX, TLorentzVector rf, TLorentzVector seconddaughter)
{
//in the case of normal to the piplus+piminus plane angles in the b1 decay the daughter = piplus, parent = omega, InverseOfX = b1, rf = gammap, seconddaughter = piminus
// boost all to rf
  TLorentzVector daughter_rf = daughter;
  TLorentzVector seconddaughter_rf = seconddaughter;
  TLorentzVector parent_rf = parent;
  TLorentzVector InverseOfX_rf = InverseOfX;
  TVector3 rfboost = rf.BoostVector();
  InverseOfX_rf.Boost(-1.0*rfboost);
  parent_rf.Boost(-1.0*rfboost);
  daughter_rf.Boost(-1.0*rfboost);
  seconddaughter_rf.Boost(-1.0*rfboost);

  //Boost daughters and parent to X
  TLorentzVector daughter_x = daughter_rf;
  TLorentzVector seconddaughter_x = seconddaughter_rf;
  TLorentzVector parent_x = parent_rf;
  TVector3 xboost = InverseOfX_rf.BoostVector();
  daughter_x.Boost(-1.0*xboost);
  seconddaughter_x.Boost(-1.0*xboost);
  parent_x.Boost(-1.0*xboost);
  
  //boost daughters to parent
  TLorentzVector daughter_parent = daughter_x;
  TLorentzVector seconddaughter_parent = seconddaughter_x;
  TVector3 parentboost = parent_x.BoostVector();
  daughter_parent.Boost(-1.0*parentboost);
  seconddaughter_parent.Boost(-1.0*parentboost);
  TVector3 daughter_parentv = daughter_parent.Vect();
  TVector3 seconddaughter_parentv = seconddaughter_parent.Vect();
  
  //get the unit vectors in space
  TVector3 daughter_parentunit = (daughter_parentv).Unit();
  TVector3 seconddaughter_parentunit = (seconddaughter_parentv).Unit();
  TVector3 normal_parentunit = daughter_parentunit.Cross(seconddaughter_parentunit);
  TVector3 InverseOfX_rfunit = (InverseOfX_rf.Vect()).Unit();

  TVector3 z = (parent_x.Vect()).Unit();
  TVector3 y = ((InverseOfX_rfunit).Cross(z)).Unit();
  TVector3 x = (y.Cross(z)).Unit();
  
  TVector3 Angles(normal_parentunit.Dot(x),normal_parentunit.Dot(y),normal_parentunit.Dot(z));

  double theta = Angles.Theta();
  double phi = Angles.Phi();

  vector <double> thetaphi{theta, phi};
    
  return thetaphi;
  
}
vector <double> getomegapiAngles(double polAngle, TLorentzVector daughter, TLorentzVector parent, TLorentzVector InverseOfX, TLorentzVector rf)
{
//in the case of omega angles in b1 decay the daughter = omega, parent = b1, InverseOfX = beam, rf = gammap
// boost all to rf
  TLorentzVector daughter_rf = daughter;
  TLorentzVector parent_rf = parent;
  TLorentzVector InverseOfX_rf = InverseOfX;
  TVector3 rfboost = rf.BoostVector();
  InverseOfX_rf.Boost(-1.0*rfboost);
  parent_rf.Boost(-1.0*rfboost);
  daughter_rf.Boost(-1.0*rfboost);

  //boost daughter to parent
  TVector3 parentboost = parent_rf.BoostVector();
  TLorentzVector daughter_parent = daughter_rf;
  daughter_parent.Boost(-1.0*parentboost);
  TVector3 daughter_parentv = daughter_parent.Vect();
  
  //get the unit vectors in space
  TVector3 daughter_parentunit = (daughter_parentv).Unit();
  TVector3 InverseOfX_rfunit = (InverseOfX_rf.Vect()).Unit();

  TVector3 z = (parent_rf.Vect()).Unit();
  TVector3 y = ((InverseOfX_rfunit).Cross(z)).Unit();
  TVector3 x = (y.Cross(z)).Unit();
  
  TVector3 Angles(daughter_parentunit.Dot(x),daughter_parentunit.Dot(y),daughter_parentunit.Dot(z));

  double theta = Angles.Theta();
  double phi = Angles.Phi();

  TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
  double Phi = atan2(y.Dot(eps), InverseOfX.Vect().Unit().Dot(eps.Cross(y)));

  vector <double> thetaphiPhi{theta, phi, Phi};
    
  return thetaphiPhi;
  
}
