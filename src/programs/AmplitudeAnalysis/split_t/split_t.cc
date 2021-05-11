#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <utility>
#include <iostream>

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"

#include "TLorentzVector.h"

#include "TH1F.h"
using namespace std;

#define DEFTREENAME "kin"

void Usage()
{
  cout << "Usage:\n  split_t <infile> <outputBase> <lowT> <highT> <nBins> [OPTIONS]\n\n";
  cout << "  Options: \n";
  cout << "   -M [maxEvents] : Limit total number of events\n";
  cout << "   -T [treeName]  : Overwrites the default ROOT tree name (\"kin\") in output and/or input files\n";
  cout << "                    To specify input and output names delimit with \':\' ex. -T inKin:outKin\n";
  cout << "   -t [treeName]  : Update existing files with new tree, instead of overwriting.\n";
  cout << "   -e             : Use a logarithmic binning (cannot have lowT = 0).\n";
  exit(1);
}


pair <string,string> GetTreeNames(char* treeArg)
{
  pair <string,string> treeNames(DEFTREENAME,"");
  string treeArgStr(treeArg);
  size_t delimPos=treeArgStr.find(':',1);

  if (delimPos != string::npos){
    treeNames.first=treeArgStr.substr(0,delimPos);
    treeNames.second=treeArgStr.substr(delimPos+1);
  }else
    treeNames.second=treeArgStr;

  return treeNames;
}


int main( int argc, char* argv[] ){
  
  unsigned int maxEvents = 4294967000; //close to 4byte int range
  
  //string treeName( "kin" );
  pair <string,string> treeNames(DEFTREENAME,DEFTREENAME);

  bool recreate=true;

  bool exponential=false;

  if( argc < 6 ) Usage();
  
  string outBase( argv[2] );
  
  double lowT = atof( argv[3] );
  double highT = atof( argv[4] );
  int numBins = atoi( argv[5] );
  
  // A somewhat convoluted way to allow tree name specification
  // via "-t [name]" in the arg. list after the standard args
  if( argc > 6 ) {
    for(int i=6; i<argc ; ++i){
      string arg=argv[i];
      if (arg == "-t"){
	if ((i+1 == argc) || (argv[i+1][0] == '-')) Usage();
	else{
	  treeNames = GetTreeNames(argv[++i]);
	  recreate=false;
	}
      }else if (arg == "-T"){
	if ((i+1 == argc) || (argv[i+1][0] == '-')) Usage();
	else{
	  treeNames = GetTreeNames(argv[++i]);
	  recreate=true;
	}
      }else if (arg == "-M"){
	if ((i+1 == argc) || (argv[i+1][0] == '-')) Usage();
	else{
	  maxEvents = atoi( argv[++i] );
	}
      }else if (arg == "-e"){
	if (lowT == 0) Usage();
	else exponential = true;
      }
      else Usage();
    }
  }

  
  vector< string > dataReaderArgs;
  dataReaderArgs.push_back( argv[1] );
  dataReaderArgs.push_back( treeNames.first );
  
  // open reader
  ROOTDataReader in( dataReaderArgs );
  
  enum { kMaxBins = 1000 };
  assert( numBins < kMaxBins );
  
  double step;
  if (exponential)
    step = ( TMath::Log10(highT) - TMath::Log10(lowT) ) / numBins;
  else
    step = ( highT - lowT ) / numBins;

  int events[numBins];
  double Tsum[numBins];
  double TsumSq[numBins+1];

  ROOTDataWriter* outFile[kMaxBins];
  
  for( int i = 0; i < numBins; ++i ){
    events[i]=0;
    Tsum[i]=0;
    TsumSq[i]=0;

    ostringstream outName;
    outName << outBase << "_" << i << ".root";
    outFile[i] = new ROOTDataWriter( outName.str(),
				     treeNames.second.c_str(), 
				     recreate, in.hasWeight());
  }
  
  unsigned int eventCount = 0;
  
  Kinematics* event;
  while( ( event = in.getEvent() ) != NULL && eventCount++ < maxEvents ){
    
    vector< TLorentzVector > fs = event->particleList();
    
    // the second entry in this list is the recoil
    TLorentzVector Target(0,0,0,0.938272046);
    TLorentzVector Recoil(fs[1]);

    double t = -1 * (Recoil - Target).M2();
    
    int bin;
    if (exponential)
      bin = static_cast< int >( floor( ( TMath::Log10(t) - TMath::Log10(lowT) ) / step ) );
    else
      bin = static_cast< int >( floor( ( t - lowT ) / step ) );
    if( ( bin < numBins ) && ( bin >= 0 ) ){
      Tsum[bin]+=t;
      TsumSq[bin]+=t*t;
      events[bin]++;
      outFile[bin]->writeEvent( *event );
      delete event;
    }
  }
  
  for( int i = 0; i < numBins; ++i ){
    double mean = Tsum[i]/events[i];
    double RMS = sqrt(TsumSq[i]/events[i] - mean*mean);
    printf("bin %2i  %10i events  mean T %.3f  RMS %.3f\n",
	   i,events[i],mean,RMS);
    delete outFile[i];
  }
  
  return 0;
}
