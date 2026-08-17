#include <map>
#include <string>
#include <stdexcept>

#include "AMPTOOLS_DATAIO/TwoPsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Histogram2D.h"
#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector2.h"
#include "TMath.h"

// Particle numbering (matches kin->particle(i)):
//   [0] = beam, [1] = Lambda, [2] = Ks, [3] = pi+
// Display labels come in through the constructor (m_partName); object names
// (the 2nd arg to the Histogram ctors) stay plain-ASCII identifiers and never
// depend on m_partName.

// ---- parse the polarization info once, at construction ----
// Walk every reaction and record how to build the beam polarization vector:
//   7-arg Zlm  -> fixed pol angle (degrees) taken from config arg 6
//   5-arg Zlm  -> pol info lives in the beam photon px/py in the tree
void TwoPsPlotGenerator::cacheArgs(){

   for( auto reaction : cfgInfo()->reactionList() ){
      const std::string reactionName = reaction->reactionName();
      const std::vector<std::string>& factor =
         cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0);
      int nargs = factor.size();

      PolInfo info;
      if( nargs == 7 ){
         info.polInPhotonP4 = false;
         info.polAngleDeg   = stod( factor.at(6) );   // real config pol angle
      } else if( nargs == 5 ){
         info.polInPhotonP4 = true;                    // pol carried in photon P4
      } else {
         throw std::invalid_argument(
            "[TwoPsPlotGenerator] Zlm amplitude has " + std::to_string(nargs) +
            " arguments (expected 5 or 7)" );
      }
      m_polInfo.emplace( reactionName, info );
   }
}

TwoPsPlotGenerator::TwoPsPlotGenerator( const FitResults& results, Option opt,
                                        const vector<string>& partName ) :
   PlotGenerator( results, opt ),      // <-- forward the option
   m_partName( partName )
{
   cacheArgs();   // populate m_polInfo before the event loop

   // ---- label helpers: TLatex display titles built from m_partName ----
   auto M2 = [&]( int i, int j ){
      return "M(" + m_partName[i] + " " + m_partName[j] + ")";
   };
   auto M1 = [&]( int i ){
      return "M(" + m_partName[i] + ")";
   };

   // resonance = particle 2 + particle 3  (K*(892) -> Ks pi+)
   const string resLabel = m_partName[2] + " " + m_partName[3];

   // item 4: use pi (not a literal 3.2) for phi ranges. TMath::Pi() dodges any
   // PI macro that might arrive from an AmpTools header.
   const double piR = TMath::Pi();

   // ---- invariant masses / Dalitz ----  (item 4: ";"-style axis titles)
   bookHistogram( kMass12, new Histogram1D( 63, 0.5, 2.0,   "hMass12", ";" + M2(1,2) + " [GeV/c^{2}]" ) );
   bookHistogram( kMass13, new Histogram1D( 63, 0.0, 2.0,   "hMass13", ";" + M2(1,3) + " [GeV/c^{2}]" ) );
   bookHistogram( kMass23, new Histogram1D( 63, 0.634, 2.203, "hMass23", ";" + M2(2,3) + " [GeV/c^{2}]" ) );
   bookHistogram( kMass1,  new Histogram1D( 60, 1.08, 1.20,   "hMass1",  ";" + M1(1) + " [GeV/c^{2}]" ) );
   bookHistogram( kMass2,  new Histogram1D( 60, 0.35, 0.65,   "hMass2",  ";" + M1(2) + " [GeV/c^{2}]" ) );
   bookHistogram( kMass3,  new Histogram1D( 60, 0.0, 0.2,   "hMass3",  ";" + M1(3) + " [GeV/c^{2}]" ) );
   bookHistogram( kDalitz, new Histogram2D( 80, 0.0, 25.0, 80, 0.0, 9.0, "Dalitz",
                     ";M^{2}(" + m_partName[1] + " " + m_partName[2] + ") [GeV^{2}/c^{4}]"
                     ";M^{2}(" + m_partName[2] + " " + m_partName[3] + ") [GeV^{2}/c^{4}]" ) );

   // ------------------------
   // ------ lab frame -------
   // ------------------------
   bookHistogram( kCosTheta_lab,     new Histogram1D( 36, -1.1, 1.1, "CosTheta_lab", ";cos#theta_{" + resLabel + "} (lab)" ) );
   bookHistogram( kPhi_lab,          new Histogram1D( 36, -piR, piR, "Phi_lab", ";#phi_{" + resLabel + "} (lab) [rad]" ) );
   bookHistogram( kCosTheta_m23_lab, new Histogram2D( 63, 0.634, 2.203, 36, -1.0, 1.0, "CosTheta_m23_lab", ";" + M2(2,3) + ";cos#theta_{lab}" ) );
   bookHistogram( kPhi_m23_lab,      new Histogram2D( 63, 0.634, 2.203, 36, -piR, piR, "Phi_m23_lab", ";" + M2(2,3) + ";#phi_{lab} [rad]" ) );

   // ------------------------
   // ---- helicity frame ----
   // ------------------------
   // ---- Polarization angle Phi (event-global, not specific to P2+P3) ----
   bookHistogram( kBigPhi_hel,       new Histogram1D( 36, -piR, piR, "BigPhi_hel", ";#Phi (helicity frame) [rad]" ) );
   // item 5: same angle with the pol. angle forced to zero -- diagnostic.
   bookHistogram( kBigPhiOffset_hel, new Histogram1D( 36, -piR, piR, "BigPhiOffset_hel", ";#Phi_{offset} (pol. angle = 0) [rad]" ) );

   // ---- Mass 2 (Ks) histograms ----
   bookHistogram( kCosThetaM2_hel,       new Histogram1D( 36, -1.1, 1.1, "CosThetaM2_hel", ";cos#theta_{" + m_partName[2] + "} (helicity frame)" ) );
   bookHistogram( kPhiM2_hel,            new Histogram1D( 36, -piR, piR, "PhiM2_hel", ";#phi_{" + m_partName[2] + "} (helicity frame) [rad]" ) );
   bookHistogram( kPhiMinusBigPhiM2_hel, new Histogram1D( 36, -piR, piR, "PhiMinusBigPhiM2_hel", ";#phi_{" + m_partName[2] + "} - #Phi (helicity frame) [rad]" ) );
   bookHistogram( kCosThetaM2_m23_hel,   new Histogram2D( 63, 0.634, 2.203, 36, -1.0, 1.0, "CosThetaM2_m23_hel", ";" + M2(2,3) + ";cos#theta_{" + m_partName[2] + "} (helicity)" ) );
   bookHistogram( kCosThetaM2_Phi_hel,   new Histogram2D( 36, -piR, piR, 36, -1.0, 1.0, "CosThetaM2_Phi_hel", ";#phi_{" + m_partName[2] + "} [rad];cos#theta_{" + m_partName[2] + "} (helicity)" ) );

   // ---- Mass 3 (pi+) histograms ----
   bookHistogram( kCosThetaM3_hel,       new Histogram1D( 36, -1.1, 1.1, "CosThetaM3_hel", ";cos#theta_{" + m_partName[3] + "} (helicity frame)" ) );
   bookHistogram( kPhiM3_hel,            new Histogram1D( 36, -piR, piR, "PhiM3_hel", ";#phi_{" + m_partName[3] + "} (helicity frame) [rad]" ) );
   bookHistogram( kPhiMinusBigPhiM3_hel, new Histogram1D( 36, -piR, piR, "PhiMinusBigPhiM3_hel", ";#phi_{" + m_partName[3] + "} - #Phi (helicity frame) [rad]" ) );
   bookHistogram( kCosThetaM3_m23_hel,   new Histogram2D( 63, 0.634, 2.203, 36, -1.0, 1.0, "CosThetaM3_m23_hel", ";" + M2(2,3) + ";cos#theta_{" + m_partName[3] + "} (helicity)" ) );
   bookHistogram( kCosThetaM3_Phi_hel,   new Histogram2D( 36, -piR, piR, 36, -1.0, 1.0, "CosThetaM3_Phi_hel", ";#phi_{" + m_partName[3] + "} [rad];cos#theta_{" + m_partName[3] + "} (helicity)" ) );

   // ---- item 5: phi vs Phi correlation, real and offset ----
   bookHistogram( kPhiM2_BigPhi_hel,       new Histogram2D( 36, -piR, piR, 36, -piR, piR, "PhiM2_BigPhi_hel",       ";#phi_{" + m_partName[2] + "} [rad];#Phi [rad]" ) );
   bookHistogram( kPhiM2_BigPhiOffset_hel, new Histogram2D( 36, -piR, piR, 36, -piR, piR, "PhiM2_BigPhiOffset_hel", ";#phi_{" + m_partName[2] + "} [rad];#Phi_{offset} [rad]" ) );

}

   std::string TwoPsPlotGenerator::numToString( Hist_index kHistName ){
    switch(kHistName){
      case kMass12: return "hMass12";
      case kMass13: return "hMass13";
      case kMass23: return "hMass23";
      case kMass1:  return "hMass1";
      case kMass2:  return "hMass2";
      case kMass3:  return "hMass3";
      case kDalitz: return "Dalitz";
      case kCosTheta_lab: return "CosTheta_lab";
      case kPhi_lab: return "Phi_lab";
      case kCosTheta_m23_lab: return "CosTheta_m23_lab";
      case kPhi_m23_lab: return "Phi_m23_lab";
      case kBigPhi_hel: return "BigPhi_hel";
      case kBigPhiOffset_hel: return "BigPhiOffset_hel";
      case kCosThetaM2_hel: return "CosThetaM2_hel";
      case kPhiM2_hel: return "PhiM2_hel";
      case kPhiMinusBigPhiM2_hel: return "PhiMinusBigPhiM2_hel";
      case kCosThetaM2_m23_hel: return "CosThetaM2_m23_hel";
      case kCosThetaM2_Phi_hel: return "CosThetaM2_Phi_hel";
      case kCosThetaM3_hel: return "CosThetaM3_hel";
      case kPhiM3_hel: return "PhiM3_hel";
      case kPhiMinusBigPhiM3_hel: return "PhiMinusBigPhiM3_hel";
      case kCosThetaM3_m23_hel: return "CosThetaM3_m23_hel";
      case kCosThetaM3_Phi_hel: return "CosThetaM3_Phi_hel";
      case kPhiM2_BigPhi_hel: return "PhiM2_BigPhi_hel";
      case kPhiM2_BigPhiOffset_hel: return "PhiM2_BigPhiOffset_hel";
      default: return "Unknown";
      }
   }

void TwoPsPlotGenerator::projectEvent( Kinematics* kin ){
   // backwards-compatible with older AmpTools (v0.10.x and prior); cannot
   // recover the lab polarization plane when multiple orientations are used.
   projectEvent( kin, "" );
}

void TwoPsPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

   // item 1 + 2: polarization comes from the cached args, not a per-event parse.
   // Empty reactionName (old single-orientation path) falls back to the first.
   const PolInfo& pol = m_polInfo.count( reactionName )
                           ? m_polInfo.at( reactionName )
                           : m_polInfo.begin()->second;

   TVector3 eps;          // beam polarization vector
   TLorentzVector P0;     // beam 4-vector used to build the angles

   if( pol.polInPhotonP4 ){
      // 5-arg Zlm: pol info stored in the beam photon px/py of the tree
      P0.SetPxPyPzE( 0., 0., kin->particle(0).E(), kin->particle(0).E() );
      eps.SetXYZ( kin->particle(0).Px(), kin->particle(0).Py(), 0. );
   } else {
      // 7-arg Zlm: fixed polarization angle (degrees) from the config -- the
      // REAL value now flows through (item 1; no more hardcoded zero).
      double a = pol.polAngleDeg * TMath::DegToRad();
      eps.SetXYZ( cos(a), sin(a), 0.0 );
      P0 = kin->particle(0);
   }

   TLorentzVector P1 = kin->particle(1);   // Lambda
   TLorentzVector P2 = kin->particle(2);   // Ks
   TLorentzVector P3 = kin->particle(3);   // pi+

   fillHistogram( kMass12, (P1+P2).M() );
   fillHistogram( kMass13, (P1+P3).M() );
   fillHistogram( kMass23, (P2+P3).M() );
   fillHistogram( kMass1,  (P1).M() );
   fillHistogram( kMass2,  (P2).M() );
   fillHistogram( kMass3,  (P3).M() );
   fillHistogram( kDalitz, (P1+P2).M2(), (P2+P3).M2() );

   TLorentzVector resonance = P2 + P3;
   TLorentzRotation resRestBoost( -resonance.BoostVector() );

  //  TLorentzVector beam_res   = resRestBoost * P0;   // (used only if you switch to GJ z-axis below.  that would require a whole refactor of the helicity frame code below.)
   TLorentzVector recoil_res = resRestBoost * P1;
   TLorentzVector p2_res     = resRestBoost * P2;
   TLorentzVector p3_res     = resRestBoost * P3;

   // ------------------------
   // ---- Lab frame ----
   // ------------------------
   GDouble locCosTheta_lab = resonance.CosTheta();
   GDouble locPhi_lab      = resonance.Phi();

   fillHistogram( kCosTheta_lab,     locCosTheta_lab );
   fillHistogram( kPhi_lab,          locPhi_lab );
   fillHistogram( kCosTheta_m23_lab, resonance.M(), locCosTheta_lab );
   fillHistogram( kPhi_m23_lab,      resonance.M(), locPhi_lab );

   // ------------------------
   // ---- helicity frame ----
   // ------------------------
   // z-axis opposite the recoil in the resonance rest frame (helicity convention).
   // For Gottfried-Jackson, use  z = beam_res.Vect().Unit()  instead.

   TVector3 z = -1. * recoil_res.Vect().Unit();
   TVector3 y = ( P0.Vect().Unit().Cross( -P1.Vect().Unit() ) ).Unit();
   TVector3 x = y.Cross(z);

   // Polarization angle Phi (real) and its zero-angle reference.
   // If kBigPhi_hel and kBigPhiOffset_hel differ, the config pol angle is being
   // applied; if they're identical, it isn't (or is genuinely zero).
   GDouble Phi = atan2( y.Dot(eps), P0.Vect().Unit().Dot( eps.Cross(y) ) );
   TVector3 epsOffset( 1., 0., 0. );   // polarization vector at angle = 0
   GDouble PhiOffset = atan2( y.Dot(epsOffset), P0.Vect().Unit().Dot( epsOffset.Cross(y) ) );

   fillHistogram( kBigPhi_hel,       Phi );
   fillHistogram( kBigPhiOffset_hel, PhiOffset );

   // Calculate and fill histograms for Mass 2 (P2 = Ks, the analyzer)
   TVector3 anglesM2( p2_res.Vect().Dot(x),
                      p2_res.Vect().Dot(y),
                      p2_res.Vect().Dot(z) );

   Double_t cosThetaM2 = anglesM2.CosTheta();
   Double_t phiM2      = anglesM2.Phi();

   fillHistogram( kCosThetaM2_hel,         cosThetaM2 );
   fillHistogram( kPhiM2_hel,              phiM2 );
   fillHistogram( kPhiMinusBigPhiM2_hel,   TVector2::Phi_mpi_pi( phiM2 - Phi ) );
   fillHistogram( kCosThetaM2_m23_hel,     (P2+P3).M(), cosThetaM2 );
   fillHistogram( kCosThetaM2_Phi_hel,     phiM2, cosThetaM2 );
   fillHistogram( kPhiM2_BigPhi_hel,       phiM2, Phi );
   fillHistogram( kPhiM2_BigPhiOffset_hel, phiM2, PhiOffset );

   // Calculate and fill histograms for Mass 3 (P3 = pi+, mirror cross-check)
   TVector3 anglesM3( p3_res.Vect().Dot(x),
                      p3_res.Vect().Dot(y),
                      p3_res.Vect().Dot(z) );

   Double_t cosThetaM3 = anglesM3.CosTheta();
   Double_t phiM3      = anglesM3.Phi();

   fillHistogram( kCosThetaM3_hel,        cosThetaM3 );
   fillHistogram( kPhiM3_hel,             phiM3 );
   fillHistogram( kPhiMinusBigPhiM3_hel,  TVector2::Phi_mpi_pi( phiM3 - Phi ) );
   fillHistogram( kCosThetaM3_m23_hel,    (P2+P3).M(), cosThetaM3 );
   fillHistogram( kCosThetaM3_Phi_hel,    phiM3, cosThetaM3 );
}
