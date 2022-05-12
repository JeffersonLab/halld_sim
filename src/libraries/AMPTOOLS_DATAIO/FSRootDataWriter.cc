#include <vector>
#include <cassert>

#include "AMPTOOLS_DATAIO/FSRootDataWriter.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"

FSRootDataWriter::FSRootDataWriter( unsigned int numParticles, const string& outFile ){
   assert(numParticles < 50);

   TH1::AddDirectory( kFALSE );
   gSystem->Load( "libTree" );

   m_outFile = new TFile( outFile.c_str(), "recreate" );
   m_outTree = new TTree( "nt", "nt" );

   m_numParticles = numParticles;

   for (unsigned int i = 0; i < m_numParticles; i++){
      TString sI("");  sI += (i+1);
      TString sEnPi = "EnP"+sI;
      TString sPxPi = "PxP"+sI;
      TString sPyPi = "PyP"+sI;
      TString sPzPi = "PzP"+sI;
      m_outTree->Branch( sEnPi, &m_EnP[i], sEnPi+"/D" );
      m_outTree->Branch( sPxPi, &m_PxP[i], sPxPi+"/D" );
      m_outTree->Branch( sPyPi, &m_PyP[i], sPyPi+"/D" );
      m_outTree->Branch( sPzPi, &m_PzP[i], sPzPi+"/D" );
   }

   m_outTree->Branch( "s12", &m_s12, "s12/D" );
   m_outTree->Branch( "s23", &m_s23, "s23/D" );

   m_outTree->Branch( "weight", &m_weight, "weight/D" );

   m_eventCounter = 0;

}


FSRootDataWriter::~FSRootDataWriter(){

   m_outFile->cd();
   m_outTree->Write();
   m_outFile->Close();

}


void
FSRootDataWriter::writeEvent( const Kinematics& kin ){

   vector< TLorentzVector > particleList = kin.particleList();

   for (unsigned int i = 0; i < m_numParticles; i++){
      m_EnP[i] = particleList[i].E();
      m_PxP[i] = particleList[i].Px();
      m_PyP[i] = particleList[i].Py();
      m_PzP[i] = particleList[i].Pz();
   }

   m_s12 = (particleList[0]+particleList[1]).M2();
   m_s23 = (particleList[1]+particleList[2]).M2();

   m_weight = kin.weight();

   m_outTree->Fill();

   m_eventCounter++;

}
