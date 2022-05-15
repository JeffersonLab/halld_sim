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

   m_outTree->Branch( "EnPB", &m_EnPB, "EnPB/D" );
   m_outTree->Branch( "PxPB", &m_PxPB, "PxPB/D" );
   m_outTree->Branch( "PyPB", &m_PyPB, "PyPB/D" );
   m_outTree->Branch( "PzPB", &m_PzPB, "PzPB/D" );
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

   m_EnPB = particleList[0].E();
   m_PxPB = particleList[0].Px();
   m_PyPB = particleList[0].Py();
   m_PzPB = particleList[0].Pz();

   for (unsigned int i = 0; i < m_numParticles; i++){
      m_EnP[i] = particleList[i+1].E();
      m_PxP[i] = particleList[i+1].Px();
      m_PyP[i] = particleList[i+1].Py();
      m_PzP[i] = particleList[i+1].Pz();
   }

   m_s12 = (particleList[0]+particleList[1]).M2();
   m_s23 = (particleList[1]+particleList[2]).M2();

   m_weight = kin.weight();

   m_outTree->Fill();

   m_eventCounter++;

}
