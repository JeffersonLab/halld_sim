
#include <vector>
#include <cassert>
#include <iostream>

#include "TLorentzVector.h"

#include "AMPTOOLS_DATAIO/ROOTDataReaderFilter.h"
#include "IUAmpTools/Kinematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include <set>

using namespace std;

ROOTDataReaderFilter::ROOTDataReaderFilter( const vector< string >& args ):
   UserDataReader< ROOTDataReaderFilter >( args ),
   m_eventCounter( 0 ),
   m_useWeight( false )
{
    
   // arguments must come in triplets (variable, min, max) with a maximum number designated by maxselects
   //    any additional multiple of 3 is possible, i.e. 1,4,7,10..., maxvar where there is an additional +1 is needed to shift to the first actual argument

   // USAGE IN AMPTOOLS CFG FILES:
   // The following amptools config line will load the data from mydata.root and perform
   //    a *SELECTION* on the Mass branch between 1.0 and 1.1. An exclamation point before
   //    the variable name will perform a *CUT*, in this case removing everything below 8.2.
   //    The -999 unphysical lower bound is there to essentially indicate no lower bound.
   // 
   // data reaction ROOTDataReaderFilter mydata.root Mass 1.0 1.1 !EBeam -999 8.2

   std::set<int> availbleNargs={1};
   for(int i=1; i<maxselects+1; ++i){ 
       availbleNargs.insert(1+3*i); }
   nargs=(int)args.size();
   nselects=(nargs-1)/3;
   assert( nselects <= maxselects );
   assert( availbleNargs.find(nargs) != availbleNargs.end() );

   string tmp;
   for (int i=0; i<nselects; ++i){
        // if the variable does not start with ! then the defined region is a selection 
        tmp=args[1+3*i];
        s_min[i] = stof(args[2+3*i]);
        s_max[i] = stof(args[3+3*i]);
        if (tmp[0]!='!'){
            tmp=args[1+3*i];
            b_isSelection[i]=true; 
        }
        // if the variable starts with ! then the defined region is a cut 
        else{
            tmp=tmp.substr(1, tmp.size());
            b_isSelection[i]=false;
        }
        s_var[i] = tmp;
        if (mapVars.find(tmp) == mapVars.end())
            mapVars[tmp]=0;
   }

   TH1::AddDirectory( kFALSE );

   //this way of opening files works with URLs of the form
   // root://xrootdserver/path/to/myfile.root
   m_inFile = TFile::Open( args[0].c_str() );

  // The name is always kin for us anyways
   m_inTree = dynamic_cast<TTree*>( m_inFile->Get( "kin" ) );

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

   m_RangeSpecified=false;
   if ( nselects >0 ){
      m_RangeSpecified=true;

      for (auto &p: mapVars){
          //cout << p.first << ' ' << p.second << endl;
          m_inTree->SetBranchAddress( p.first.c_str(), &p.second );
      }

      cout << "*********************************************" << endl;
      cout << "NVARS: " << nselects << endl;
      cout << "Total events: " <<  m_inTree->GetEntries() << endl;
      
      m_numEvents = 0;
      m_weightIntegral = 0;
      for (int i=0; i<nselects; ++i){
          if (b_isSelection[i])
            cout << "Selecting " << s_min[i] << " < " << s_var[i] << " < " << s_max[i] << endl;
          else
            cout << "Cutting " << s_min[i] << " < " << s_var[i] << " < " << s_max[i] << endl;
      }
      
      while( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
         m_inTree->GetEntry( m_eventCounter++ );
         assert( m_nPart < Kinematics::kMaxParticles );
      
         vector< TLorentzVector > particleList;
         particleList.push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );
         for( int i = 0; i < m_nPart; ++i ){
            particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
         }
      
         selection=true;
         for (int i=0; i<nselects; ++i){
            if (b_isSelection[i]){
                // For a selection: we cut an event if the values are less than the minimum or greater than the maximum
                if ( ( (mapVars[s_var[i]]<s_min[i]) || (mapVars[s_var[i]]>s_max[i]) ) ){
                    selection=false;            
                    break;
                } 
            }
            else{
                // For a cut: we cut an event if the values are in between the maximum and minimum
                if ( ( (mapVars[s_var[i]]>s_min[i]) && (mapVars[s_var[i]]<s_max[i]) ) ){
                    selection=false;            
                    break;
                } 
            }
         }
      
        if (selection){
            m_numEvents++;
            m_weightIntegral+=m_weight;
        }
      }
      cout << "[" << args[0] << "] Number of events kept    = " << m_numEvents << endl;
      cout << "[" << args[0] << "] Weight Integral    = " << m_weightIntegral << endl;
      cout << "*********************************************" << endl;
   }
}

ROOTDataReaderFilter::~ROOTDataReaderFilter()
{
   if( m_inFile != NULL ) m_inFile->Close();
}

void ROOTDataReaderFilter::resetSource()
{

   cout << "Resetting source " << m_inTree->GetName() 
      << " in " << m_inFile->GetName() << endl;

   // this will cause the read to start back at event 0
   m_eventCounter = 0;
}

   Kinematics*
ROOTDataReaderFilter::getEvent()
{
   if (m_RangeSpecified == false){
      if( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
         m_inTree->GetEntry( m_eventCounter++ );
         assert( m_nPart < Kinematics::kMaxParticles );

         vector< TLorentzVector > particleList;
         particleList.push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );

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
         particleList.push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );

         for( int i = 0; i < m_nPart; ++i ){
            particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
         }

         selection=true;
         for (int i=0; i<nselects; ++i){
            if (b_isSelection[i]){
                // For a selection: we cut an event if the values are less than the minimum or greater than the maximum
                if ( ( (mapVars[s_var[i]]<s_min[i]) || (mapVars[s_var[i]]>s_max[i]) ) ){
                    selection=false;            
                    break;
                } 
            }
            else{
                // For a cut: we cut an event if the values are in between the maximum and minimum
                if ( ( (mapVars[s_var[i]]>s_min[i]) && (mapVars[s_var[i]]<s_max[i]) ) ){
                    selection=false;            
                    break;
                } 
            }
         }
         if (!selection) continue;

         return new Kinematics( particleList, m_useWeight ? m_weight : 1.0 ); 
      }
      return NULL;
   }
}

unsigned int ROOTDataReaderFilter::numEvents() const
{	
   return m_numEvents;
}