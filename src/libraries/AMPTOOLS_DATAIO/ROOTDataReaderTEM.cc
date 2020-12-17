
#include <vector>
#include <cassert>
#include <iostream>

#include "TLorentzVector.h"

#include "ROOTDataReaderTEM.h"
#include "IUAmpTools/Kinematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

ROOTDataReaderTEM::ROOTDataReaderTEM( const vector< string >& args ):
   UserDataReader< ROOTDataReaderTEM >( args ),
   m_eventCounter( 0 ),
   m_useWeight( false )
{
   assert( args.size() == 8 || args.size() == 7 ); //TEM cuts with and without special tree name

   TH1::AddDirectory( kFALSE );

   //this way of opening files works with URLs of the form
   // root://xrootdserver/path/to/myfile.root
   m_inFile = TFile::Open( args[0].c_str() );

   // default to tree name of "kin" if none is provided
   if( args.size() == 8 ){

      m_inTree = dynamic_cast<TTree*>( m_inFile->Get( args[7].c_str() ) );
   }
   else{

      m_inTree = dynamic_cast<TTree*>( m_inFile->Get( "kin" ) );
   }

   m_numEvents = m_inTree->GetEntries();

   m_inTree->SetBranchAddress( "NumFinalState", &m_nPart );
   m_inTree->SetBranchAddress( "E_FinalState", m_e );
   m_inTree->SetBranchAddress( "Px_FinalState", m_px );
   m_inTree->SetBranchAddress( "Py_FinalState", m_py );
   m_inTree->SetBranchAddress( "Pz_FinalState", m_pz );
   m_inTree->SetBranchAddress( "E_Beam", &m_eBeam );
   m_inTree->SetBranchAddress( "Px_Beam", &m_pxBeam );
   m_inTree->SetBranchAddress( "Py_Beam", &m_pyBeam );
   m_inTree->SetBranchAddress( "Pz_Beam", &m_pzBeam );

   if(m_inTree->GetBranch("Weight") != NULL){

     m_useWeight = true;
     m_inTree->SetBranchAddress( "Weight", &m_weight );
   }
   else{
     
     m_useWeight=false;
   }

   m_RangeSpecified = false;
   if( args.size() == 8 || args.size() == 7){
      // Set t range
      m_tMin = atof(args[1].c_str());
      m_tMax = atof(args[2].c_str());
      m_EMin = atof(args[3].c_str());
      m_EMax = atof(args[4].c_str());
      m_MMin = atof(args[5].c_str());
      m_MMax = atof(args[6].c_str());
      m_RangeSpecified = true;

      m_numEvents = 0;
      cout << "*********************************************" << endl;
      cout << "ROOT Data reader  -t range specified [" << m_tMin << "," << m_tMax << ")" << endl;
      cout << "ROOT Data reader Beam E range specified [" << m_EMin << "," << m_EMax << ")" << endl;
      cout << "ROOT Data reader  Inv. Mass range specified [" << m_MMin << "," << m_MMax << ")" << endl;
      cout << "Total events: " <<  m_inTree->GetEntries() << endl;

      while( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){

         m_inTree->GetEntry( m_eventCounter++ );
         assert( m_nPart < Kinematics::kMaxParticles );

         vector< TLorentzVector > particleList;
         TLorentzVector finalstate;
	 TLorentzVector recoil;

         particleList.
            push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );

         for( int i = 0; i < m_nPart; ++i ){

            particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
	    
	    if (i > 0 && i < 5) finalstate += TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] );
	    if (i == 0 || i == 5) recoil += TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] );
         }

         // Calculate -t and check if it is in range
         // Use the reconstructed proton
         TLorentzVector target = TLorentzVector(0.0,0.0,0.0,0.938272);
         double tMag = fabs((target-recoil).M2());
         double EMag = TLorentzVector(m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ).E();

	 double MMag = finalstate.M();

         if (m_tMin <= tMag && tMag < m_tMax && m_EMin <= EMag && EMag < m_EMax && m_MMin <= MMag && MMag < m_MMax){
            m_numEvents++;
         }
      }
      cout << "Number of events kept    = " << m_numEvents << endl;
      cout << "*********************************************" << endl;
   }   
}

ROOTDataReaderTEM::~ROOTDataReaderTEM()
{
   if( m_inFile != NULL ) m_inFile->Close();
}

void ROOTDataReaderTEM::resetSource()
{

   cout << "Resetting source " << m_inTree->GetName() 
      << " in " << m_inFile->GetName() << endl;

   // this will cause the read to start back at event 0
   m_eventCounter = 0;
}

   Kinematics*
ROOTDataReaderTEM::getEvent()
{


   if (m_RangeSpecified == false){

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

         return new Kinematics( particleList, m_useWeight ? m_weight : 1.0 );
      }
      else return NULL;

      } 
      else{

         while( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){

            m_inTree->GetEntry( m_eventCounter++ );
            assert( m_nPart < Kinematics::kMaxParticles );

            vector< TLorentzVector > particleList;
	    TLorentzVector finalstate;
	    TLorentzVector recoil;

            particleList.
               push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );

            for( int i = 0; i < m_nPart; ++i ){

               particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
	       if (i > 0 && i < 5) finalstate += TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] );
	       if (i == 0 || i == 5) recoil += TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] );
            }

            // Calculate -t and check if it is in range
            // Use the reconstructed proton
            TLorentzVector target = TLorentzVector(0.0,0.0,0.0,0.938272);
            double tMag = fabs((target-recoil).M2());
	    double EMag = TLorentzVector(m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ).E();

	    double MMag = finalstate.M();

         if (m_tMin <= tMag && tMag < m_tMax && m_EMin <= EMag && EMag < m_EMax && m_MMin <= MMag && MMag < m_MMax){
               return new Kinematics( particleList, m_useWeight ? m_weight : 1.0 ); 
            }

         }
         return NULL;
      }

      return NULL;
   }

   unsigned int ROOTDataReaderTEM::numEvents() const
   {	
      return m_numEvents;
   }
