// decay_evtgen.cc
// Description
// Sean Dobbs, sdobbs@fsu.edu (2019)

#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
using namespace std;

#include "HDDM/hddm_s.hpp"

#include "EVTGEN_WRAPPER/EvtGenInterface.h"


string INPUT_FILE = "";
string OUTPUT_FILE = "";
//string USER_DECAY = "userDecay.dec";

bool PROCESS_ALL_EVENTS = true;
int NUM_EVENTS_TO_PROCESS = -1;
//bool GEN_SCHANNEL = false;

void ParseCommandLineArguments(int narg,char *argv[], EvtGenInterface *evt_gen_interface);
void Usage(void);



//-------------------------------
// main
//-------------------------------ƒ
int main(int narg, char *argv[])
{
	EvtGenInterface *evt_gen_interface = new EvtGenInterface();

	ParseCommandLineArguments(narg,argv,evt_gen_interface);

	if (INPUT_FILE == "") {
	  cerr << "No input file!" << endl;
	}

	// Open input file
	ifstream *infile = new ifstream(INPUT_FILE);
	if (! infile->is_open()) {
	  cerr << "Unable to open file \"" << INPUT_FILE << "\" for reading."
				<< endl;
	  exit(-2);
	}
	cout << "Opening Input File:  " << INPUT_FILE << " ..." << endl;
	hddm_s::istream *instream = new hddm_s::istream(*infile);

	// Open output file
	ofstream *outfile = new ofstream(OUTPUT_FILE.c_str());
	if (! outfile->is_open()) {
	  cerr << "Unable to open output file \"" << OUTPUT_FILE
				<< "\" for writing." << endl;
	  exit(-3);
	}
	cout << "Opening Output File:  " << OUTPUT_FILE << " ..." << endl;
	hddm_s::ostream *outstream = new hddm_s::ostream(*outfile);

	evt_gen_interface->InitEvtGen();
	
	int event_count = 1;
	hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
	while(*instream >> *hddmevent) {
	  // next line commented out, an unused variable
	  //		int num_particles = -1;   // number of particles in the event
		if( (event_count++%1000) == 0) {
			cout << "Processed " << event_count << " events ..." << endl;
		}

		evt_gen_interface->Decay(hddmevent);

	   	*outstream << *hddmevent;  // save event

		// see if we should stop processing
		if(!PROCESS_ALL_EVENTS) {
		  if(event_count >= NUM_EVENTS_TO_PROCESS)
		    break;
		}
	}

	// cleanup
	delete instream;
	delete outstream;

	return 0;
}

//-------------------------------
// ConvertStringInt
//-------------------------------
bool ConvertStringInt(string &the_str, int &out_int) 
{
	istringstream the_istream(the_str);
	the_istream >> out_int;
	if(the_istream.fail())
		return false;

	return true;
}

//-------------------------------
// ParseCommandLineArguments
//-------------------------------
void ParseCommandLineArguments(int narg,char *argv[], EvtGenInterface *evt_gen_interface)
{
  string num_events_str;
  size_t seperator_index;
  string argstr;
  int pdgtype1, pdgtype2;
  string substr1, substr2;

   if (narg < 2) {
      Usage();
      exit(0);
   }

   for(int i=1; i < narg; i++) {
      if (argv[i][0]=='-') {
         char *ptr = &argv[i][1];
         switch(*ptr) {
            case 'n':
	      PROCESS_ALL_EVENTS = false;
	      num_events_str = &ptr[1];
	      NUM_EVENTS_TO_PROCESS = std::stoi(num_events_str);
              break;
            case 'o':
              OUTPUT_FILE = &ptr[1];
              break;
            case 'u':
              argstr = &ptr[1];
              evt_gen_interface->SetUserDecay(argstr);
              break;
            case 'S':
              evt_gen_interface->SetGenSChannel(true);
              break;
            case 'X':
              // assume input of the form -XN_M
              // where N and M are both PDG ID numbers
              argstr = &ptr[1];
              seperator_index = argstr.find("_");
              if(seperator_index == string::npos) {
              	cerr << " Invalid -X format: " << argstr << endl;
              } else {
              	// let's actually try parsing
              	substr1 = argstr.substr(0,seperator_index);
              	if(!ConvertStringInt(substr1,pdgtype1)) {
              		cerr << " Invalid -X format: " << argstr << " - bad type " << substr1 << endl;
              	} else {
              		substr2 = argstr.substr(seperator_index+1);
              		if(!ConvertStringInt(substr2,pdgtype2)) {
              	 		cerr << " Invalid -X format: " << argstr << " - bad type " << substr2 << endl;
              	 	} else {
              	 		cerr << "setting up conversion of type " << pdgtype1 << " to " << pdgtype2 << endl;
              	 		evt_gen_interface->AddPDGConversionMap(pdgtype1, pdgtype2);
              	 	}
              	 }
              }
              break;
            default:
              cerr << "Unknown option \"" << argv[i] << "\"" << endl;
              Usage();
              exit(-1);
         }
      }
      else {
         INPUT_FILE = argv[i];
      }
   }
   
   if(OUTPUT_FILE == "") {
	   // Determine output filename from input filename
	   OUTPUT_FILE = INPUT_FILE;
	   size_t pos = OUTPUT_FILE.find_last_of(".");
	   if (pos != string::npos) OUTPUT_FILE.erase(pos);
	   OUTPUT_FILE += "_decayed.hddm";
   }
}

//-------------------------------
// Usage
//-------------------------------
void Usage(void)
{
  cout << endl;
  cout << "Usage:" << endl;
  cout << "       decay_evtgen [options] file.hddm" << endl;
  cout << endl;
  cout << "Decay thrown particles via EvtGen" << endl;
  cout << "The particle types and decays are defined in $EVTGENDIR/evt.pdl " << endl
       << "  and $EVTGENDIR/DECAY.DEC respectively." << endl;
  cout << "The location of these files can be override by the environment variables " << endl
       << "  EVTGEN_PARTICLE_DEFINITIONS and EVTGEN_DECAY_FILE respectively." << endl;
  cout << endl << "For more documentation on EvtGen, go to https://evtgen.hepforge.org/" << endl;
  cout << endl;
  cout << " options:" << endl;
  cout << endl;
  cout << "  -nNumEvents               "
               "number of events to process (default: all)" << endl;
  cout << "  -o\"output_file_name\"    "
               "set the file name used for output (default: append \"_decayed\")" << endl;
  cout << "  -u\"user_decay_file_name\"    "
               "set the file name of the user decay file (default: userDecay.dec)" << endl;
  cout << "  -XM_N                         "
      "translate particles of PDG ID M to those of ID N.  This option can be specified multiple times" << endl;
  cout << "  -h                        "
               "print this usage statement." << endl;
  cout << endl;
}

