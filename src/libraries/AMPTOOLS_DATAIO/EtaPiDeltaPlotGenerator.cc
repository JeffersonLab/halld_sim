#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/EtaPiDeltaPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

EtaPiDeltaPlotGenerator::EtaPiDeltaPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
    // calls to bookHistogram go here
    
    // bookHistogram( k2PiMass, new Histogram1D( 200, 0., 2.0, "M2pi", "Invariant Mass of #pi^{+} #pi^{-}") );
    bookHistogram( kEtaPiMass, new Histogram1D( 80, 0.4, 3.0, "EtaPiMass", "Invariant Mass of #eta #pi^{-}") );
    bookHistogram( kDeltaPPMass, new Histogram1D( 60, 0.8, 2.0, "DeltaPPMass", "Invariant Mass of p#pi^{+}") );
    
    bookHistogram( kEtaCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );
    bookHistogram( kPhi, new Histogram1D( 50, -180, 180, "Phi", "#Phi" ) );
    bookHistogram( kt, new Histogram1D( 100, 0, 2.00, "t", "-t" ) );
}

void
EtaPiDeltaPlotGenerator::projectEvent( Kinematics* kin ){
    
    TLorentzVector beam   = kin->particle( 0 );
    TLorentzVector protonP4 = kin->particle( 1 );//proton
    TLorentzVector p1 = kin->particle( 2 ); //Eta
    TLorentzVector p2 = kin->particle( 3 ); //Pi-
    TLorentzVector p3 = kin->particle( 4 ); //Pi+
    
    TLorentzVector recoil=protonP4+p3;
    
    TLorentzVector resonance = p1 + p2;
    TLorentzRotation resonanceBoost( -resonance.BoostVector() );
    
    TLorentzVector recoil_res = resonanceBoost * recoil;
    TLorentzVector p1_res = resonanceBoost * p1;
    
    TLorentzVector locCoMP4=recoil + resonance;
    TVector3 boostCoM=-(locCoMP4.Vect())*(1.0/locCoMP4.E());

    TLorentzVector locBeamP4_CM=beam;
    TLorentzVector locEtaP4_CM=p1;
    TLorentzVector locPiMinusP4_CM=p2;
    TLorentzVector locEtaPiMinusP4_CM=resonance;
    TLorentzVector locDeltaPPP4_CM=recoil;
    
    locDeltaPPP4_CM.Boost(boostCoM);
    locEtaPiMinusP4_CM.Boost(boostCoM);
    locBeamP4_CM.Boost(boostCoM);
    locEtaP4_CM.Boost(boostCoM);
    locPiMinusP4_CM.Boost(boostCoM);
    
    //GJ Boost
    TVector3 boostGJ=-(locEtaPiMinusP4_CM.Vect())*(1.0/locEtaPiMinusP4_CM.E());
    
    //Define GJ frame vectors
    TLorentzVector locEtaPiMinus_GJ=locEtaPiMinusP4_CM;
    TLorentzVector locEtaP4_GJ=locEtaP4_CM;
    TLorentzVector locPiminusP4_GJ=locEtaP4_CM;
    TLorentzVector locBeamP4GJ=locBeamP4_CM;
    
    //Boost in GJ
    locEtaPiMinus_GJ.Boost(boostGJ);
    locBeamP4GJ.Boost(boostGJ);
    locEtaP4_GJ.Boost(boostGJ);
    locPiminusP4_GJ.Boost(boostGJ);
    
    TVector3 z_GJ;
    z_GJ.SetXYZ(locBeamP4GJ.X(),locBeamP4GJ.Y(),locBeamP4GJ.Z());//z GJ
    TVector3 z_hat_GJ=z_GJ.Unit();
    TVector3 y_GJ=locBeamP4_CM.Vect().Cross(locEtaPiMinusP4_CM.Vect());
    TVector3 y_hat_GJ=y_GJ.Unit();
    TVector3 x_hat_GJ=y_hat_GJ.Cross(z_hat_GJ);//x hat GJ
    
    TVector3 v(locEtaP4_GJ.Vect()*x_hat_GJ, locEtaP4_GJ.Vect()*y_hat_GJ,
               locEtaP4_GJ.Vect()*z_hat_GJ);
    double cosTheta = v.CosTheta();
    // comment out next line, unused variable
    //    double theta = v.Theta();
    double phi = v.Phi()*180./TMath::Pi();
    
    // normal to the production plane
    //TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

    // compute invariant t
    //GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());
    TLorentzVector TargetP4;
    TargetP4.SetPxPyPzE(0,0,0,0.938272);
    GDouble t=(recoil-TargetP4).Mag2();
    // calls to fillHistogram go here
    
    fillHistogram( kEtaPiMass, ( resonance ).M() );
    
    fillHistogram( kDeltaPPMass, ( recoil ).M() );
    fillHistogram( kEtaCosTheta, cosTheta );
    fillHistogram( kPhi, phi );
    fillHistogram( kt, -t );      // fill with -t to make positive
}
