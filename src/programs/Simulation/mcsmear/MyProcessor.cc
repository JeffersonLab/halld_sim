// Author: David Lawrence  Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.cc
//

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <JANA/Compatibility/JGeometryManager.h>


using namespace std;

#include <strings.h>

#include "MyProcessor.h"
#include "hddm_s_merger.h"
#include "DANA/DEvent.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <JANA/Compatibility/JGeometryXML.h>
#include <TRACKING/DMCThrown.h>
#include <DRandom2.h>
#include <FDCSmearer.h>

extern char *OUTFILENAME;
extern std::map<hddm_s::istream*,double> files2merge;
extern std::map<hddm_s::istream*,hddm_s::streamposition> start2merge;
extern std::map<hddm_s::istream*,int> skip2merge;

static pthread_mutex_t output_file_mutex;
static pthread_t output_file_mutex_last_owner;
static pthread_mutex_t input_file_mutex;
static pthread_t input_file_mutex_last_owner;

#include <JANA/Calibrations/JCalibration.h>
//static JCalibration *jcalib=NULL;
static bool locCheckCCDBContext = true;

static thread_local Smear *smearer(0);
static pthread_mutex_t smearer_mutex;
static pthread_t smearer_mutex_last_owner;

//-----------
// PrintCCDBWarning
//-----------
static void PrintCCDBWarning(string context)
{
	jout << endl;
	jout << "===============================================================================" << endl;
	jout << "                            !!!!! WARNING !!!!!" << endl;
	jout << "You have either not specified a CCDB variation, or specified a variation" << endl;
	jout << "that appears inconsistent with processing simulated data." << endl;
	jout << "Be sure that this is what you want to do!" << endl;
	jout << endl;
	jout << "  JANA_CALIB_CONTEXT = " << context << endl;
	jout << endl;
	jout << "For a more detailed list of CCDB variations used for simulations" << endl;
	jout << "see the following wiki page:" << endl;
	jout << endl;
	jout << "   https://halldweb.jlab.org/wiki/index.php/Simulations#Simulation_Conditions" << endl;
	jout << "===============================================================================" << endl;
	jout << endl;
}



void mcsmear_thread_HUP_sighandler(int sig)
{
   jerr<<" Caught HUP signal for thread 0x"<<hex<<pthread_self()<<dec<<" thread exiting..."<<endl;

   // We use output_file_mutex_owner to keep track (sort of)
   // of which thread has the mutex locked. This mutex is only
   // locked at the end of the evnt method. Once the lock is
   // obtained, this value is set to hold the id of the thread
   // that locked it. It may help in debugging to know the last
   // known thread to have locked the mutex when the signal
   // handler was called
   jerr<<endl;
   jerr<<" Last thread to lock output file mutex: 0x"<<hex<<pthread_self()<<dec<<endl;
   jerr<<" Attempting to unlock mutex to avoid deadlock." <<endl;
   jerr<<" However, the output file is likely corrupt at" <<endl;
   jerr<<" this point and the process should be restarted ..." <<endl;
   jerr<<endl;
   pthread_mutex_unlock(&output_file_mutex);
   pthread_exit(NULL);
}


//------------------------------------------------------------------
// Init   -Open output file 
//------------------------------------------------------------------
void MyProcessor::Init(void)
{
   // open HDDM file
   ofs = new ofstream(OUTFILENAME);
   if (!ofs->is_open()){
      cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
      exit(-1);
   }
   fout = new hddm_s::ostream(*ofs);
   Nevents_written = 0;

   HDDM_USE_COMPRESSION = 2;
   auto app = GetApplication();
   app->SetDefaultParameter("HDDM:USE_COMPRESSION", HDDM_USE_COMPRESSION,
                          "Turn on/off compression of the output HDDM stream."
                          " \"0\"=no compression, \"1\"=bz2 compression, \"2\"=z compression (default)");
   HDDM_USE_INTEGRITY_CHECKS = true;
   app->SetDefaultParameter("HDDM:USE_INTEGRITY_CHECKS",
                                HDDM_USE_INTEGRITY_CHECKS,
                          "Turn on/off automatic integrity checking on the"
                          " output HDDM stream."
                          " Set to \"0\" to turn off (it's on by default)");

   // enable on-the-fly bzip2 compression on output stream
   if (HDDM_USE_COMPRESSION == 0) {
      jout << " HDDM compression disabled" << std::endl;
   } else if (HDDM_USE_COMPRESSION == 1) {
      jout << " Enabling bz2 compression of output HDDM file stream" 
           << std::endl;
      fout->setCompression(hddm_s::k_bz2_compression);
   } else {
      jout << " Enabling z compression of output HDDM file stream (default)" 
           << std::endl;
      fout->setCompression(hddm_s::k_z_compression);
   }

   // enable a CRC data integrity check at the end of each event record
   if (HDDM_USE_INTEGRITY_CHECKS) {
      jout << " Enabling CRC data integrity check in output HDDM file stream"
           << std::endl;
      fout->setIntegrityChecks(hddm_s::k_crc32_integrity);
   }
   else {
      jout << " HDDM integrity checks disabled" << std::endl;
   }

   // We set the mutex type to "ERRORCHECK" so that if the
   // signal handler is called, we can unlock the mutex
   // safely whether we have it locked or not.
   pthread_mutexattr_t attr;
   pthread_mutexattr_init(&attr);
   pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK);
   pthread_mutex_init(&output_file_mutex, NULL);
   pthread_mutex_init(&input_file_mutex, NULL);
   pthread_mutex_init(&smearer_mutex, NULL);
   
   // pthreads does not provide an "invalid" value for 
   // a pthread_t that we can initialize with. Furthermore,
   // the pthread_t may be a simple as an integer or as
   // a complicated structure. Hence, to make this portable
   // we clear it with bzero.
   bzero(&output_file_mutex_last_owner, sizeof(pthread_t));
   bzero(&input_file_mutex_last_owner, sizeof(pthread_t));
   bzero(&smearer_mutex_last_owner, sizeof(pthread_t));
   
   // By default, the empirical fdc DOCA-dependent efficiency
   // is applied to fdc wire hits in FDCSmearer. This can be
   // disabled by this commandline option.
   app->SetDefaultParameter("FDC:USE_EFFVSDOCA", fdc_config_t::FDC_EFFVSDOCA,
                          "Turn on/off fdc doca-dependent wire hit efficiency,"
                          " \"0\"=no efficiency correction,"
                          " \"1\"=standard efficiency correction (default)");

   return; //NOERROR;
}

void MyProcessor::BeginRun(const std::shared_ptr<const JEvent>& event)
{
    // Generally, simulations should be generated and analyzed with a non-default
    // set of calibrations, since the calibrations needed for simulations are
    // different than those needed for data.  
    // Conventionally, the CCDB variations needed for simulations start with the 
    // string "mc".  To guard against accidentally not setting the variation correctly
    // we check to see if the variation is set and if it contains the string "mc".
    // Note that for now, we only print a warning and do not exit immediately.
    // It might be advisable to apply some tougher love.
   auto locRunNumber = event->GetRunNumber();
    if(locCheckCCDBContext) {
        // only do this once per brun record
        locCheckCCDBContext = true;

        // load the CCDB context
        auto app =  event->GetJApplication();
        auto run_number = event->GetRunNumber();
        JCalibration* jcalib = app->GetService<JCalibrationManager>()->GetJCalibration(run_number);
    
        string context = jcalib->GetContext();
      
        jout << "checking context = " << context << endl;

        // Really we should parse the context string, but since "mc" shouldn't show up
        // outside of the context, we just search the whole string.
        // Also make sure that the variation is being set
        if( (context.find("variation") == string::npos) || (context.find("mc") == string::npos) ) {
            PrintCCDBWarning(context);
        }

        std::map<string, float> parms;
        jcalib->Get("TOF/tof_parms", parms);
        hddm_s_merger::set_ftof_min_delta_t_ns(parms.at("TOF_TWO_HIT_RESOL"));
        jcalib->Get("FDC/fdc_parms", parms);
        hddm_s_merger::set_fdc_wires_min_delta_t_ns(parms.at("FDC_TWO_HIT_RESOL"));
        jcalib->Get("START_COUNTER/start_parms", parms);
        hddm_s_merger::set_stc_min_delta_t_ns(parms.at("START_TWO_HIT_RESOL"));
        jcalib->Get("BCAL/bcal_parms", parms);
        hddm_s_merger::set_bcal_min_delta_t_ns(parms.at("BCAL_TWO_HIT_RESOL"));
        jcalib->Get("FCAL/fcal_parms", parms);
        hddm_s_merger::set_fcal_min_delta_t_ns(parms.at("FCAL_TWO_HIT_RESOL"));
    }
   

	// load configuration parameters for all the detectors
    pthread_mutex_lock(&smearer_mutex);
    smearer_mutex_last_owner = pthread_self();

	if(smearer != NULL)
		delete smearer;
	smearer = new Smear(config, event, config->DETECTORS_TO_LOAD);

#ifdef HAVE_RCDB
	// Pull configuration parameters from RCDB
	bool haveRCDBConfigFile = false;
	if(!config->SKIP_READING_RCDB) {
	  haveRCDBConfigFile = config->ParseRCDBConfigFile(locRunNumber);
	}
	if(haveRCDBConfigFile) {

	        const double fadc250_period_ns(4.);
	        const double fadc125_period_ns(8.);
		
		// Default parameters

		int cdc_npeak     =  1;
		double cdc_ie     =  200.;
		double cdc_pg     =  4;

		int fdc_nhits     =  100;
		int fdc_npeak     =  1;
		double fdc_width  =  80.;
		double fdc_ie     =  16;
		double fdc_pg     =  4;

		int stc_npeak     =  3;
		int stc_nhits     =  8;
		double stc_width  =  80;
		double stc_nsa    =  20;
		double stc_nsb    =  5;
             
		int bcal_npeak     =  1;
		int bcal_nhits     =  8;
		double bcal_width  =  100;
		double bcal_nsa    =  26;
		double bcal_nsb    =  1;

		int tof_npeak     =  3;
		int tof_nhits     =  64;
		double tof_width  =  80;
		double tof_nsa    =  10;
		double tof_nsb    =  1;

		int fcal_npeak     =  3;
		double fcal_nsa    =  15;
		double fcal_nsb    =  1;

		int ps_npeak   =  3;
		double ps_nsa  =  10;
		double ps_nsb  =  3;

		int psc_npeak     =  3;
		int psc_nhits     =  8;
		double psc_nsa    =  6;
		double psc_nsb    =  3;
		double psc_width  =  100;

		int tagm_npeak     =  3;
		int tagm_nhits     =  8;

		double tagm_width  =  60;
		double tagm_nsa    =  6;
		double tagm_nsb    =  3;

		int tpol_npeak  =  1;


		// hits merging / truncation parameters for the CDC
		if(config->readout["CDC"].size() > 0){
		  cdc_npeak   =  config->readout["CDC"].at("NPEAK");
		  cdc_ie      =  config->readout["CDC"].at("IE");
		  cdc_pg      =  config->readout["CDC"].at("PG");
		}

		double cdc_gate = (cdc_ie + cdc_pg) * fadc125_period_ns;

		hddm_s_merger::set_cdc_max_hits(cdc_npeak);
		hddm_s_merger::set_cdc_integration_window_ns(cdc_gate);
		


		// hits merging / truncation parameters for the FDC
		if(config->readout["FDC"].size() > 0){
		  fdc_nhits  =  config->readout["FDC"].at("NHITS");
		  fdc_npeak  =  config->readout["FDC"].at("NPEAK");
		  fdc_width  =  config->readout["FDC"].at("WIDTH");
		  fdc_ie     =  config->readout["FDC"].at("IE");
		  fdc_pg     =  config->readout["FDC"].at("PG");
		}

		double fdc_gate = (fdc_ie + fdc_pg) * fadc125_period_ns;

		hddm_s_merger::set_fdc_wires_max_hits(fdc_nhits);
		hddm_s_merger::set_fdc_wires_min_delta_t_ns(fdc_width + 5.);
		hddm_s_merger::set_fdc_strips_max_hits(fdc_npeak);
		hddm_s_merger::set_fdc_strips_integration_window_ns(fdc_gate);



		// hits merging / truncation parameters for the STC
		if(config->readout["ST"].size() > 0){
		  stc_npeak = config->readout["ST"].at("NPEAK");
		  stc_nhits = config->readout["ST"].at("NHITS");
		  stc_width = config->readout["ST"].at("WIDTH");
		  stc_nsa = config->readout["ST"].at("NSA");
		  stc_nsb = config->readout["ST"].at("NSB");
		}
		
		double stc_gate = (stc_nsa + stc_nsb) * fadc250_period_ns;

		hddm_s_merger::set_stc_adc_max_hits(stc_npeak);
		hddm_s_merger::set_stc_tdc_max_hits(stc_nhits);
		hddm_s_merger::set_stc_min_delta_t_ns(stc_width + 5.);
		hddm_s_merger::set_stc_integration_window_ns(stc_gate);
		

		// hits merging / truncation parameters for the BCAL
		if(config->readout["BCAL"].size() > 0){
		  bcal_npeak  =  config->readout["BCAL"].at("NPEAK");
		  bcal_nhits  =  config->readout["BCAL"].at("NHITS");
		  bcal_width = config->readout["BCAL"].at("WIDTH");
		  bcal_nsa = config->readout["BCAL"].at("NSA");
		  bcal_nsb = config->readout["BCAL"].at("NSB");
		}

		double bcal_gate = (bcal_nsa + bcal_nsb) * fadc250_period_ns;

		hddm_s_merger::set_bcal_adc_max_hits(bcal_npeak);
		hddm_s_merger::set_bcal_tdc_max_hits(bcal_nhits);
		hddm_s_merger::set_bcal_min_delta_t_ns(bcal_width + 5.);
		hddm_s_merger::set_bcal_integration_window_ns(bcal_gate);

	       
		
		// hits merging / truncation parameters for the TOF
		if(config->readout["TOF"].size() > 0){
		  tof_npeak = config->readout["TOF"].at("NPEAK");
		  tof_nhits = config->readout["TOF"].at("NHITS");
		  tof_width = config->readout["TOF"].at("WIDTH");
		  tof_nsa = config->readout["TOF"].at("NSA");
		  tof_nsb = config->readout["TOF"].at("NSB");		  
		}

		double tof_gate = (tof_nsa + tof_nsb) * fadc250_period_ns;

		hddm_s_merger::set_ftof_adc_max_hits(tof_npeak);
		hddm_s_merger::set_ftof_tdc_max_hits(tof_nhits);
		hddm_s_merger::set_ftof_min_delta_t_ns(tof_width + 5.);
		hddm_s_merger::set_ftof_integration_window_ns(tof_gate);
		

		// hits merging / truncation parameters for the FCAL
		if(config->readout["FCAL"].size() > 0){
		  fcal_npeak = config->readout["FCAL"].at("NPEAK");
		  fcal_nsb = config->readout["FCAL"].at("NSB");
		  
		}

		double fcal_gate = (fcal_nsa + fcal_nsb) * fadc250_period_ns;

		hddm_s_merger::set_fcal_max_hits(fcal_npeak);
		hddm_s_merger::set_fcal_integration_window_ns(fcal_gate);


		// hits merging / truncation parameters for the ECAL
		int ecal_npeak     =  3;
		double ecal_nsa    =  15;
		double ecal_nsb    =  1;
		
		if(config->readout["ECAL"].size() > 0){
		  ecal_npeak = config->readout["ECAL"].at("NPEAK");
		  ecal_nsb = config->readout["ECAL"].at("NSB");		  
		}
		
		double ecal_gate = (ecal_nsa + ecal_nsb) * fadc250_period_ns;
		
		hddm_s_merger::set_ecal_max_hits(ecal_npeak);
		hddm_s_merger::set_ecal_integration_window_ns(ecal_gate);
				
		
		// hits merging / truncation parameters for the CCAL
		int ccal_npeak     =  3;
		double ccal_nsa    =  15;
		double ccal_nsb    =  1;
		
		if(config->readout["CCAL"].size() > 0){
		  ccal_npeak = config->readout["CCAL"].at("NPEAK");
		  ccal_nsb = config->readout["CCAL"].at("NSB");		  
		}

		double ccal_gate = (ccal_nsa + ccal_nsb) * fadc250_period_ns;
		
		hddm_s_merger::set_ccal_max_hits(ccal_npeak);
		hddm_s_merger::set_ccal_integration_window_ns(ccal_gate);

		


		// hits merging / truncation parameters for the PS
		if(config->readout["PS"].size() > 0){
		  ps_npeak = config->readout["PS"].at("NPEAK");
		  ps_nsa = config->readout["PS"].at("NSA");
		  ps_nsb = config->readout["PS"].at("NSB");
		}

		if(config->readout["PSC"].size() > 0){
		  psc_npeak = config->readout["PSC"].at("NPEAK");
		  psc_nhits = config->readout["PSC"].at("NHITS");
		  psc_nsa = config->readout["PSC"].at("NSA");
		  psc_nsb = config->readout["PSC"].at("NSB");
		  psc_width = config->readout["PSC"].at("WIDTH");
		}

		double ps_gate  = (ps_nsa + ps_nsb) * fadc250_period_ns;
		double psc_gate = (psc_nsa + psc_nsb) * fadc250_period_ns;

		hddm_s_merger::set_ps_max_hits(ps_npeak);
		hddm_s_merger::set_ps_integration_window_ns(ps_gate);

		hddm_s_merger::set_psc_adc_max_hits(psc_npeak);
		hddm_s_merger::set_psc_tdc_max_hits(psc_nhits);
		hddm_s_merger::set_psc_min_delta_t_ns(psc_width + 5.);
		hddm_s_merger::set_psc_integration_window_ns(psc_gate);
		

		// hits merging / truncation parameters for the TAGM/TAGH
		if(config->readout["TAGM"].size() > 0){
		  tagm_npeak = config->readout["TAGM"].at("NPEAK");
		  tagm_nhits = config->readout["TAGM"].at("NHITS");
		  tagm_width = config->readout["TAGM"].at("WIDTH");
		  tagm_nsa = config->readout["TAGM"].at("NSA");
		  tagm_nsb = config->readout["TAGM"].at("NSB");
		}

		double tagm_gate = (tagm_nsa + tagm_nsb) * fadc250_period_ns;
		
		hddm_s_merger::set_tag_adc_max_hits(tagm_npeak);
		hddm_s_merger::set_tag_tdc_max_hits(tagm_nhits);
		hddm_s_merger::set_tag_min_delta_t_ns(tagm_width + 5.);
		hddm_s_merger::set_tag_integration_window_ns(tagm_gate);
		

		// hits merging / truncation parameters for the TPOL		
		if(config->readout["TPOL"].size() > 0){
		  tpol_npeak = config->readout["TPOL"].at("NPEAK");
		}
		
		hddm_s_merger::set_tpol_max_hits(tpol_npeak);


	}

#endif  // HAVE_RCDB

    // fast forward any merger input files over skipped events
    std::map<hddm_s::istream*,hddm_s::streamposition>::iterator iter;
    for (iter = start2merge.begin(); iter != start2merge.end(); ++iter) {
        hddm_s::HDDM record2;

        iter->first->setPosition(start2merge.at(iter->first));
        iter->first->skip(skip2merge[iter->first]);
          
		if (!(*iter->first >> record2)) {
			std::cerr << "Trying to merge from empty input file, "
					  << "cannot continue!" << std::endl;
			exit(-1);
		}

        skip2merge[iter->first] = 0;
    }

    hddm_s_merger::set_config_run_loaded(locRunNumber);
    pthread_mutex_unlock(&smearer_mutex);
    return; //NOERROR;
}

//------------------------------------------------------------------
// Process - Do processing for each event here
//------------------------------------------------------------------
void MyProcessor::Process(const std::shared_ptr<const JEvent>& event)
{
   JEventSource *source = event->GetJEventSource();
   DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
   if (!hddm_source) {
      cerr << " This program MUST be used with an HDDM file as input!" << endl;
      exit(-1);
   }
   hddm_s::HDDM *record = const_cast<hddm_s::HDDM*>(event->GetSingle<hddm_s::HDDM>());
   if (!record)
      return; //NOERROR;
 
   // Make sure the run-dependent config has been loaded for this run
   hddm_s::PhysicsEventList pev = record->getPhysicsEvents();
   if (pev.size() > 0) {
      if (hddm_s_merger::get_config_run_loaded() != pev(0).getRunNo()) {
         BeginRun(event);
      }
   }

   // Handle geometry records
   hddm_s::GeometryList geom = record->getGeometrys();
   if (geom.size() > 0) {
   JGeometryXML *jgeom = dynamic_cast<JGeometryXML*>(DEvent::GetDGeometry(event)->GetJGeometry());
      geom(0).setMd5smear(jgeom->GetChecksum());
      if (geom(0).getMd5smear() != geom(0).getMd5simulation()) {
         std::cerr << "Warning: simulation geometry checksum does not match"
                   << " the geometry description used by mcsmear."
                   << std::endl;
      }
   }
   
   // Smear values
   smearer->SmearEvent(record);

   // Load any external events to be merged during smearing
   std::map<hddm_s::istream*,double>::iterator iter;
   for (iter = files2merge.begin(); iter != files2merge.end(); ++ iter) {
      int count = iter->second;
      if (count != iter->second) {
         count = gDRandom.Poisson(iter->second);
      }
      for (int i=0; i < count; ++i) {
         hddm_s::HDDM record2;
         if (!(*iter->first >> record2)) {
            //pthread_mutex_lock(&input_file_mutex);
            //input_file_mutex_last_owner = pthread_self();
            iter->first->setPosition(start2merge.at(iter->first));
            if (!(*iter->first >> record2)) {
               //pthread_mutex_unlock(&input_file_mutex);
               std::cerr << "Trying to merge from empty input file, "
                         << "cannot continue!" << std::endl;
               exit(-1);
            }
            //pthread_mutex_unlock(&input_file_mutex);
         }
         
         if(config->MERGE_TAGGER_HITS == false) {
         	hddm_s_merger::set_tag_merging(false);
         }
         hddm_s_merger::set_t_shift_ns(0);
         hddm_s::RFsubsystemList RFtimes = record2.getRFsubsystems();
         hddm_s::RFsubsystemList::iterator RFiter;
         for (RFiter = RFtimes.begin(); RFiter != RFtimes.end(); ++RFiter)
            if (RFiter->getJtag() == "TAGH")
               hddm_s_merger::set_t_shift_ns(-RFiter->getTsync());
         *record += record2;
      }
   }

   // Apply DAQ truncation to hit lists
   if (config->APPLY_HITS_TRUNCATION)
      hddm_s_merger::truncate_hits(*record);

   // Write event to output file
   //pthread_mutex_lock(&output_file_mutex);
   //output_file_mutex_last_owner = pthread_self();
   *fout << *record;
   Nevents_written++;
   //pthread_mutex_unlock(&output_file_mutex);

   return; //NOERROR;
}

//------------------------------------------------------------------
// Finish   -Close output file here
//------------------------------------------------------------------
void MyProcessor::Finish(void)
{
   if (fout)
      delete fout;
   if (ofs) {
      ofs->close();
      cout << endl << "Closed HDDM file" << endl;
   }
   cout << " " << Nevents_written << " event written to " << OUTFILENAME
        << endl;
   
   return; //NOERROR;
}
