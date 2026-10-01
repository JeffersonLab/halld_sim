#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"
#include "decayAnglesTest.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

// helper function to calculate vector (secondary) decay
void calcVectorDecayAngles( VecPSAngles& outputAngles,
                            const TLorentzVector& vecVPS,       // four-vectors in the vec-ps frame
                            const TLorentzVector& vecDaught1VPS,
                            const TLorentzVector& vecDaught2VPS, 
                            const TVector3& zAxisParent ){      // z-axis of the initial decay

    // boost daughter particles to parent (vector)'s rest frame
    TVector3 vecBoost = -vecVPS.BoostVector();
    TLorentzVector vecDaught1Vec = vecDaught1VPS;
    TLorentzVector vecDaught2Vec = vecDaught2VPS;
    vecDaught1Vec.Boost( vecBoost );
    vecDaught2Vec.Boost( vecBoost );
 
    // define helicity frame of the vector particle
    TVector3 z = vecVPS.Vect().Unit();
    TVector3 y = zAxisParent.Cross( z ).Unit();  
    TVector3 x = y.Cross( z );

    // for omega->3pi decay, the decay vector is normal to the decay plane
    // otherwise, the decay vector is one of the daughters
    TVector3 decayVector = ( vecDaught2Vec.E() > 0 ? vecDaught1Vec.Vect().Cross( vecDaught2Vec.Vect() ) : vecDaught1Vec.Vect() );

    TVector3 components( decayVector.Dot( x ),
                         decayVector.Dot( y ),
                         decayVector.Dot( z ) );

    outputAngles.cosThetaH = components.CosTheta();
    outputAngles.phiH = components.Phi();

    // compute the variable lambda for omega->3pi decay
    double m0Sq = 0.0182187;    // mass of pi0, squared
    double mqSq = 0.0194798;    // mass of pi+/-, squared
    double m3piSq = vecVPS.M2();
    double lambda_max = 0.75 * TMath::Power( 1./9. * ( 5*m3piSq + 3*( m0Sq - 4*mqSq ) - 4*sqrt( m3piSq*m3piSq + 3*m3piSq*( m0Sq - mqSq ) ) ), 2 );
    outputAngles.lambda = fabs( decayVector.Mag2() ) / lambda_max;
}

// Lower vertex (recoil) angles in the Gottfried-Jackson frame
LVAngles getLVAnglesGJ( const TLorentzVector& beamLab,
                        const TLorentzVector& recoilLab,
                        const TLorentzVector& protonLab ){

    LVAngles angles;

    angles.bigPhiLV = recoilLab.Vect().Phi();

    TVector3 recoilBoost = -recoilLab.BoostVector();
    TLorentzVector protonRecoil = protonLab;
    protonRecoil.Boost( recoilBoost );

    TLorentzVector targetRecoil(0,0,0,0.9382720813);
    targetRecoil.Boost( recoilBoost );

    // TODO: document discussion of axis defns
    TVector3 z = targetRecoil.Vect().Unit();
    TVector3 y = -beamLab.Vect().Cross( recoilLab.Vect() ).Unit();
    TVector3 x = y.Cross( z );

    TVector3 components( protonRecoil.Vect().Dot( x ),
                         protonRecoil.Vect().Dot( y ),
                         protonRecoil.Vect().Dot( z ) );

    angles.cosThetaLV = components.CosTheta();
    angles.phiLV = components.Phi();

    return angles;
}

// Lower vertex (recoil) angles in the helicity frame
LVAngles getLVAnglesHelicity( const TLorentzVector& beamLab,
                              const TLorentzVector& recoilLab,
                              const TLorentzVector& protonLab){

    LVAngles angles;

    angles.bigPhiLV = recoilLab.Vect().Phi();

    // boost daughter proton to recoil frame
    TVector3 recoilBoost = -recoilLab.BoostVector();
    TLorentzVector protonRecoil = protonLab;
    protonRecoil.Boost( recoilBoost );

    TLorentzVector targetLab(0,0,0,0.9382720813);
    TLorentzVector resonanceLab = targetLab + beamLab - recoilLab;
    TLorentzVector resonanceRecoil = resonanceLab;
    resonanceRecoil.Boost( recoilBoost );

    // TODO: document discussion of axis defns
    TVector3 z = -resonanceRecoil.Vect().Unit();
    TVector3 y = -beamLab.Vect().Unit().Cross( z );
    TVector3 x = y.Cross( z );

    TVector3 components( protonRecoil.Vect().Dot( x ),
                         protonRecoil.Vect().Dot( y ),
                         protonRecoil.Vect().Dot( z ) );
    
    angles.cosThetaLV = components.CosTheta();
    angles.phiLV = components.Phi();

    return angles;
}


// VecPS two-step decay, with primary decay in the GJ frame
VecPSAngles getVecPSAnglesGJ( const TLorentzVector& beamLab,
                              const TLorentzVector& vecPSLab,
                              const TLorentzVector& vecLab,
                              const TLorentzVector& vecDaught1Lab,
                              const TLorentzVector& vecDaught2Lab ){

    VecPSAngles angles;

    // boost vector and beam to vec-ps rest frame
    TVector3 vecPSBoost = -vecPSLab.BoostVector();
    TLorentzVector vecVPS = vecLab;   // vector in the vec-ps rest frame
    vecVPS.Boost( vecPSBoost );

    TLorentzVector beamVPS = beamLab;
    beamVPS.Boost( vecPSBoost );

    TVector3 z = beamVPS.Vect().Unit();
    TVector3 y = beamLab.Vect().Cross( vecPSLab.Vect() ).Unit();
    TVector3 x = y.Cross( z );    

    TVector3 components( vecVPS.Vect().Dot( x ),
                         vecVPS.Vect().Dot( y ),
                         vecVPS.Vect().Dot( z ) );

    angles.cosThetaX = components.CosTheta();
    angles.phiX = components.Phi();

//    angles.bigPhiX = vecPSLab.Vect().Phi();
    TVector3 eps(1,0,0);
    angles.bigPhiX = atan2( y.Dot( eps ), beamLab.Vect().Unit().Dot( eps.Cross( y ) ) );

    // boost vector's daughters to vec-ps rest frame
    TLorentzVector vecDaught1VPS = vecDaught1Lab;
    TLorentzVector vecDaught2VPS = vecDaught2Lab;
    vecDaught1VPS.Boost( vecPSBoost );
    vecDaught2VPS.Boost( vecPSBoost );

    calcVectorDecayAngles( angles, vecVPS, vecDaught1VPS, vecDaught2VPS, z );

    return angles;
}

// VecPS two-step decay, with primary decay in the helicity frame
VecPSAngles getVecPSAnglesHelicity( const TLorentzVector& beamLab,
                                    const TLorentzVector& vecPSLab,
                                    const TLorentzVector& vecLab,
                                    const TLorentzVector& vecDaught1Lab,
                                    const TLorentzVector& vecDaught2Lab ){

    VecPSAngles angles;

    // boost recoil and vector to vecPS rest frame
    TLorentzVector targetLab(0,0,0,0.9382720813);
    TLorentzVector recoilLab = targetLab + beamLab - vecPSLab;
    
    TVector3 vecPSBoost = -vecPSLab.BoostVector();
    TLorentzVector vecVPS = vecLab;   // vector in the vec-ps rest frame
    vecVPS.Boost( vecPSBoost );

    TLorentzVector recoilVPS = recoilLab;
    recoilVPS.Boost( vecPSBoost );

    // recoil in the vec-ps rest frame is in the opposite direction
    // of the vec-ps in the CM rest frame, which is the direction
    // of the helicity z-axis
    // the y axis is k cross z, where k is the direction of the beam
    // in the CM rest frame. It has the same direction in the lab frame,
    // so that one is used here
    TVector3 z = -recoilVPS.Vect().Unit();
    TVector3 y = beamLab.Vect().Unit().Cross( z ).Unit();
    TVector3 x = y.Cross( z );

    TVector3 components( vecVPS.Vect().Dot( x ),
                         vecVPS.Vect().Dot( y ),
                         vecVPS.Vect().Dot( z ) );

    angles.cosThetaX = components.CosTheta();
    angles.phiX = components.Phi();

//    angles.bigPhiX = vecPSLab.Vect().Phi();
    TVector3 eps(1,0,0);
    angles.bigPhiX = atan2( y.Dot( eps ), beamLab.Vect().Unit().Dot( eps.Cross( y ) ) );

    // boost vector's daughters to vec-ps rest frame
    TLorentzVector vecDaught1VPS = vecDaught1Lab;
    TLorentzVector vecDaught2VPS = vecDaught2Lab;
    vecDaught1VPS.Boost( vecPSBoost );
    vecDaught2VPS.Boost( vecPSBoost );

    calcVectorDecayAngles( angles, vecVPS, vecDaught1VPS, vecDaught2VPS, z );

    return angles; 
}
