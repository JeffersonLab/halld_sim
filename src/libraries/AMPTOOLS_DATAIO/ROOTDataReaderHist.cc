
#include <vector>
#include <cassert>
#include <iostream>

#include "TLorentzVector.h"

#include "AMPTOOLS_DATAIO/ROOTDataReaderHist.h"
#include "IUAmpTools/Kinematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

ROOTDataReaderHist::ROOTDataReaderHist( const vector< string >& args ):
  UserDataReader< ROOTDataReaderHist >( args ),
  m_eventCounter( 0 ),
  m_useWeight( false )
{
  assert( args.size() == 3 );
  
  TH1::AddDirectory( kFALSE );
  
  //this way of opening files works with URLs of the form
  // root://xrootdserver/path/to/myfile.root
  m_inFile = TFile::Open( args[0].c_str() );
  
  // default to tree name of "kin"
  m_inTree = dynamic_cast<TTree*>( m_inFile->Get( "kin" ) );
  
  m_inTree->SetBranchAddress( "NumFinalState", &m_nPart );
  m_inTree->SetBranchAddress( "E_FinalState", m_e );
  m_inTree->SetBranchAddress( "Px_FinalState", m_px );
  m_inTree->SetBranchAddress( "Py_FinalState", m_py );
  m_inTree->SetBranchAddress( "Pz_FinalState", m_pz );
  m_inTree->SetBranchAddress( "E_Beam", &m_eBeam );
  m_inTree->SetBranchAddress( "Px_Beam", &m_pxBeam );
  m_inTree->SetBranchAddress( "Py_Beam", &m_pyBeam );
  m_inTree->SetBranchAddress( "Pz_Beam", &m_pzBeam );

  if(m_inTree->GetBranch("Weight") != NULL) {

    m_useWeight = true;
    m_inTree->SetBranchAddress( "Weight", &m_weight );
  }
  else{

    m_useWeight = false;
  }
  
  //root file name provided as argument
  m_HistFile = TFile::Open( args[1].c_str() );
  //second argument: histogram name 
  m_Hist = dynamic_cast<TH1F*>( m_HistFile->Get( args[2].c_str() ) );

}

ROOTDataReaderHist::~ROOTDataReaderHist()
{
  if( m_inFile != NULL ) m_inFile->Close();
}

void
ROOTDataReaderHist::resetSource()
{
	
  cout << "Resetting source " << m_inTree->GetName() 
       << " in " << m_inFile->GetName() << endl;
  
  // this will cause the read to start back at event 0
  m_eventCounter = 0;
}

Kinematics*
ROOTDataReaderHist::getEvent()
{
  if( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
    //  if( m_eventCounter < 10 ){
    
    m_inTree->GetEntry( m_eventCounter++ );
    assert( m_nPart < Kinematics::kMaxParticles );
    
    vector< TLorentzVector > particleList;
    
    particleList.
      push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );
    
    for( int i = 0; i < m_nPart; ++i ){
      
      particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
    }

    float mass = (particleList[2]+particleList[3]).M();
    float hist_weight = m_Hist->GetBinContent(m_Hist->FindBin(mass));
    
    return new Kinematics( particleList, m_useWeight ? m_weight*hist_weight : hist_weight );
  }
  else{
    
    return NULL;
  }
}

unsigned int
ROOTDataReaderHist::numEvents() const
{	
  return static_cast< unsigned int >( m_inTree->GetEntries() );
}
