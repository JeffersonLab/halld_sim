//The purpose of this function is to return the decay angles of the daughter particle in the helicity frame of the parent.
//For a two particle decay, the function arguments are the lab frame TLorentzVectors of: the daughter particle, the parent particle, the inverse of the X-axis and the reference Z-axis.
//For a three particle decay, the function arguments are the same with the addition of the second daughter particle and the X-axis reference.
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

#include "HelicityFrame.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"
 
vector <double>  getthetaphi(TLorentzVector particle1, TLorentzVector z, TLorentzVector InverseOfX, TLorentzVector zrf, TLorentzVector particle2, TLorentzVector xrf)
{
  TVector3 xrfboost = xrf.BoostVector();
  TLorentzVector z_rf = z;
  TLorentzVector InverseOfX_rf = InverseOfX;
  TVector3 zrfboost = zrf.BoostVector();

//   if (xrf) boost zrf, z, x to xrf
  if(xrf.E() > 0.0){
//     TVector3 xrfboost = xrf->BoostVector();
    zrf.Boost(-1.0*xrfboost);
    z_rf.Boost(-1.0*xrfboost);
    InverseOfX_rf.Boost(-1.0*xrfboost);
    zrfboost = zrf.BoostVector();//update zrf boost vector
  }
// boost x to zrf
  else{
  InverseOfX_rf.Boost(-1.0*zrfboost);
  }

  //boost z to rf
  z_rf.Boost(-1.0*zrfboost);

  //particle 1 (positively charged)
  TLorentzVector particle1_rf = particle1;
 // if (xrf) boost p1  to xrf
  if(xrf.E() > 0.0){
    particle1_rf.Boost(-1.0*xrfboost);
  }
  //boost particle 1 to zrf
  particle1_rf.Boost(-1.0*zrfboost);
  //boost particle 1 to z
  TVector3 zboost = z_rf.BoostVector();
  TLorentzVector particle1_z = particle1_rf;
  particle1_z.Boost(-1.0*zboost);
  TVector3 particle_z = particle1_z.Vect();

   
     //if (p2 exists) (negatively charged) 
     if (particle2.E() > 0.0){
  TLorentzVector particle2_rf = particle2;
 // if (xrf) boost p2  to xrf
  if(xrf.E() > 0.0){
    particle2_rf.Boost(-1.0*xrfboost);
  }
  //boost particle 2 to zrf
  particle2_rf.Boost(-1.0*zrfboost);
  //boost particle 2 to z
  //TVector3 zboost = z_rf.BoostVector();
  TLorentzVector particle2_z = particle2_rf;
  particle2_z.Boost(-1.0*zboost);
  //normal to decay plane = p2.Cross(p2)
  particle_z = (particle1_z.Vect()).Cross(particle2_z.Vect());
  }

  //get the unit vectors in space
  TVector3 particle_zunit = (particle_z).Unit();
  TVector3 z_rfunit = (z_rf.Vect()).Unit();
  TVector3 InverseOfX_rfunit = (InverseOfX_rf.Vect()).Unit();

 //calculate theta  
  double theta = particle_zunit.Angle(z_rfunit);

  //calculate phi
  TVector3 y = ((InverseOfX_rfunit).Cross(z_rfunit)).Unit();
  TVector3 x = (y.Cross(z_rfunit)).Unit();
  double phi = TMath::ATan2(particle_zunit.Dot(y), particle_zunit.Dot(x));
    
  vector <double> thetaphi{theta, phi};
    
  return thetaphi;
  
}


vector <double> getthetaphi(TLorentzVector particle1, TLorentzVector z, TLorentzVector InverseOfX, TLorentzVector zrf)
{

  TLorentzVector z_rf = z;
  TLorentzVector InverseOfX_rf = InverseOfX;
  TVector3 zrfboost = zrf.BoostVector();

// boost x to zrf
  InverseOfX_rf.Boost(-1.0*zrfboost);

  //boost z to rf
  z_rf.Boost(-1.0*zrfboost);

  //particle 1 (positively charged)
  TLorentzVector particle1_rf = particle1;

  //boost particle 1 to zrf
  particle1_rf.Boost(-1.0*zrfboost);
  //boost particle 1 to z
  TVector3 zboost = z_rf.BoostVector();
  TLorentzVector particle1_z = particle1_rf;
  particle1_z.Boost(-1.0*zboost);
  TVector3 particle_z = particle1_z.Vect();
  
  //get the unit vectors in space
  TVector3 particle_zunit = (particle_z).Unit();
  TVector3 z_rfunit = (z_rf.Vect()).Unit();
  TVector3 InverseOfX_rfunit = (InverseOfX_rf.Vect()).Unit();

 //calculate theta  
  double theta = particle_zunit.Angle(z_rfunit);

  //calculate phi
  TVector3 y = ((InverseOfX_rfunit).Cross(z_rfunit)).Unit();
  TVector3 x = (y.Cross(z_rfunit)).Unit();
  double phi = TMath::ATan2(particle_zunit.Dot(y), particle_zunit.Dot(x));

  vector <double> thetaphi{theta, phi};
    
  return thetaphi;
  
}
