#include "FCALSmearer.h"

//-----------
// fcal_config_t  (constructor)
//-----------
fcal_config_t::fcal_config_t(JEventLoop *loop, const DFCALGeometry *fcalGeom) 
{
	// default values
	FCAL_PHOT_STAT_COEF     = 0.0; // 0.05;
	FCAL_BLOCK_THRESHOLD    = 0.0; // 20.0*k_MeV;
	FCAL_TSIGMA             = 0.0; // 400 ps 
	FCAL_PED_RMS            = 0.0; // 3.0
	FCAL_MC_ESCALE          = 0.0; // 1.54
	FCAL_ADC_ASCALE         = 0.0; // 2.7E-4 
	FCAL_INTEGRAL_PEAK      = 0.0; // 5.7
	FCAL_THRESHOLD          = 0.0; // 108
	FCAL_THRESHOLD_SCALING  = 0.0; // (110/108)
	FCAL_ENERGY_WIDTH_FLOOR = 0.0; // 0.03

	// Get values from CCDB
	cout << "Get FCAL/fcal_parms parameters from CCDB..." << endl;
    map<string, double> fcalparms;
    if(loop->GetCalib("FCAL/fcal_parms", fcalparms)) { 
     	jerr << "Problem loading FCAL/fcal_parms from CCDB!" << endl;
    } else {
       	FCAL_PHOT_STAT_COEF   = fcalparms["FCAL_PHOT_STAT_COEF"]; 
       	FCAL_BLOCK_THRESHOLD  = fcalparms["FCAL_BLOCK_THRESHOLD"];
	}
		
	cout<<"get FCAL/gains from calibDB"<<endl;
    vector <double> FCAL_GAINS_TEMP;
    if(loop->GetCalib("FCAL/gains", FCAL_GAINS_TEMP)) {
    	jerr << "Problem loading FCAL/gains from CCDB!" << endl;
    } else {
    	for (unsigned int i = 0; i < FCAL_GAINS_TEMP.size(); i++) {
       		FCAL_GAINS.push_back(FCAL_GAINS_TEMP.at(i));
    	}
    }
     
   cout<<"get FCAL/pedestals from calibDB"<<endl;
   vector <double> FCAL_PEDS_TEMP;
   if(loop->GetCalib("FCAL/pedestals", FCAL_PEDS_TEMP)) {
      jerr << "Problem loading FCAL/pedestals from CCDB!" << endl;
   } else {
       for (unsigned int i = 0; i < FCAL_PEDS_TEMP.size(); i++) {
         FCAL_PEDS.push_back(FCAL_PEDS_TEMP.at(i));
      }
   }
   
   cout<<"get FCAL/MC/pedestal_rms from calibDB"<<endl;
   double FCAL_PED_RMS_TEMP;
   if(loop->GetCalib("FCAL/MC/pedestal_rms", FCAL_PED_RMS_TEMP)) {
      jerr << "Problem loading FCAL/MC/pedestal_rms from CCDB!" << endl;
   } else {
      FCAL_PED_RMS = FCAL_PED_RMS_TEMP;
   }
	
   cout<<"get FCAL/MC/integral_peak_ratio from calibDB"<<endl;
   double FCAL_INT_PEAK_TEMP;
   if(loop->GetCalib("FCAL/MC/integral_peak_ratio", FCAL_INT_PEAK_TEMP)) {
      jerr << "Problem loading FCAL/MC/integral_peak_ratio from CCDB!" << endl;
   } else {
      FCAL_INTEGRAL_PEAK = FCAL_INT_PEAK_TEMP;
   }
   
   cout<<"get FCAL/MC/threshold from calibDB"<<endl;
   double FCAL_THRESHOLD_TEMP;
   if(loop->GetCalib("FCAL/MC/threshold", FCAL_THRESHOLD_TEMP)) {
      jerr << "Problem loading FCAL/MC/threshold from CCDB!" << endl;
   } else {
      FCAL_THRESHOLD = FCAL_THRESHOLD_TEMP;
   }
   
   cout<<"get FCAL/MC/threhsold_scaling from calibDB"<<endl;
   double FCAL_THRESHOLD_SCALING_TEMP;
   if(loop->GetCalib("FCAL/MC/threshold_scaling", FCAL_THRESHOLD_SCALING_TEMP)) {
       jerr << "Problem loading FCAL/MC/threshold_scaling from CCDB!" << endl;
   } else {
       FCAL_THRESHOLD_SCALING = FCAL_THRESHOLD_SCALING_TEMP;
   }

   cout<<"get FCAL/MC/mc_escale from calibDB"<<endl;
   double FCAL_MC_ESCALE_TEMP;
   if(loop->GetCalib("FCAL/MC/mc_escale", FCAL_MC_ESCALE_TEMP)) {
        jerr << "Problem loading FCAL/MC/mc_escale from CCDB!" << endl;
   } else {
        FCAL_MC_ESCALE = FCAL_MC_ESCALE_TEMP;
   }

    cout<<"get FCAL/MC/energy_width_floor from calibDB"<<endl;
    double FCAL_ENERGY_WIDTH_FLOOR_TEMP;
    if(loop->GetCalib("FCAL/MC/energy_width_floor", FCAL_ENERGY_WIDTH_FLOOR_TEMP)) {
        jerr << "Problem loading FCAL/MC/energy_width_floor from CCDB!" << endl;
    } else {
        FCAL_ENERGY_WIDTH_FLOOR = FCAL_ENERGY_WIDTH_FLOOR_TEMP;
    }

    cout<<"get FCAL/digi_scales parameters from calibDB"<<endl;
    map<string, double> fcaldigiscales;
    if(loop->GetCalib("FCAL/MC/digi_scales", fcaldigiscales)) {
    	jerr << "Problem loading FCAL/MC/digi_scales from CCDB!" << endl;
    } else {
        FCAL_ADC_ASCALE = fcaldigiscales["FCAL_ADC_ASCALE"];
    }

    cout<<"get FCAL/mc_timing_smear parameters from calibDB"<<endl;
    map<string, double> fcalmctimingsmear;
    if(loop->GetCalib("FCAL/mc_timing_smear", fcalmctimingsmear)) {
    	jerr << "Problem loading FCAL/mc_timing_smear from CCDB!" << endl;
    } else {
        FCAL_TSIGMA = fcalmctimingsmear["FCAL_TSIGMA"];
    }



    // initialize 2D matrix of efficiencies, indexed by (row,column)
    vector< vector<double > > new_block_efficiencies(DFCALGeometry::kBlocksTall, 
						     vector<double>(DFCALGeometry::kBlocksWide));
    block_efficiencies = new_block_efficiencies;
    
    // load efficiencies from CCDB and fill 
    vector<double> raw_table;

    if(loop->GetCalib("FCAL/block_mc_efficiency", raw_table)) {
      jerr << "Problem loading FCAL/block_mc_efficiency from CCDB!" << endl;
    } else {
      for (int channel=0; channel < static_cast<int>(raw_table.size()); channel++) {
	int row = fcalGeom->row(channel);
	int col = fcalGeom->column(channel);
	block_efficiencies[row][col] = raw_table[channel];
      }
    }


    //   7/27/2020 A.S.  Exclude run-by-run determined bad channels listed in the /FCAL/block_quality table for PrimEx runs, 
    //                     to make simulation consistent with reconstruction
    //   Channels efficiency (FCAL/block_mc_efficiency) can be applied after excluding bad channels
   

    int primex_run = 0;

    if (loop->GetCalib("/PHOTON_BEAM/pair_spectrometer/experiment", primex_run))
      jerr << "Problem loading /PHOTON_BEAM/pair_spectrometer/experment/run from CCDB!" << endl;


    if(primex_run == 1){
      
      vector< double > raw_block_qualities;    // we should change this to an int?
      
      int BAD_CH = 1;
      
      
      if (loop->GetCalib("/FCAL/block_quality", raw_block_qualities))
	jout << "/FCAL/block_quality not used for this run" << endl;
      else {
	
	for (int channel=0; channel < static_cast<int>(raw_block_qualities.size()); channel++) {
	  int row = fcalGeom->row(channel);
	  int col = fcalGeom->column(channel);
	  // Exclude bad channels
	  if(raw_block_qualities[channel] == BAD_CH)
	    block_efficiencies[row][col] = -1.;  
	}
	
      }
    } 
                  
}
	
//-----------
// SmearEvent
//-----------
void FCALSmearer::SmearEvent(hddm_s::HDDM *record)
{
   /// Smear the FCAL hits using the nominal resolution of the individual blocks.
   /// The way this works is a little funny and warrants a little explanation.
   /// The information coming from hdgeant is truth information indexed by 
   /// row and column, but containing energy deposited and time. The mcsmear
   /// program will copy the truth information from the fcalTruthHit element
   /// to a new fcalHit element, smearing the values with the appropriate detector
   /// resolution.
   ///
   /// To access the "truth" values in DANA, get the DFCALHit objects using the
   /// "TRUTH" tag.

   //if (!fcalGeom)
   //   fcalGeom = new DFCALGeometry();

   hddm_s::FcalBlockList blocks = record->getFcalBlocks();
   hddm_s::FcalBlockList::iterator iter;
   for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
      iter->deleteFcalHits();
      hddm_s::FcalTruthHitList thits = iter->getFcalTruthHits();
      hddm_s::FcalTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
	int row=iter->getRow();
	int column=iter->getColumn();

         // Simulation simulates a grid of blocks for simplicity. 
         // Do not bother smearing inactive blocks. They will be
         // discarded in DEventSourceHDDM.cc while being read in
         // anyway.
         if (!fcalGeom->isBlockActive(row, column))
            continue;

	 double E = titer->getE();
	 double Ethreshold=fcal_config->FCAL_BLOCK_THRESHOLD;
	 double sigEfloor=fcal_config->FCAL_ENERGY_WIDTH_FLOOR;
	 double sigEstat=fcal_config->FCAL_PHOT_STAT_COEF;
	 double Erange = 8.0;
         
	 if (row<DFCALGeometry::kBlocksTall&&column<DFCALGeometry::kBlocksWide){
	   // correct simulation efficiencies 
	   if (config->APPLY_EFFICIENCY_CORRECTIONS
	       && !gDRandom.DecideToAcceptHit(fcal_config->GetEfficiencyCorrectionFactor(row, column))) {
	     continue;
	   } 
	 
	   // Get gain constant per block	      
	   int channelnum = fcalGeom->channel(row, column); 
	   double FCAL_gain = fcal_config->FCAL_GAINS.at(channelnum);      
	   double pedestal_rms = fcal_config->FCAL_PED_RMS;
	   double integral_peak = fcal_config->FCAL_INTEGRAL_PEAK;
	   double MeV_FADC = fcal_config->FCAL_ADC_ASCALE;
	   double pedestal = fcal_config->FCAL_PEDS.at(channelnum);
	   double threshold = fcal_config->FCAL_THRESHOLD;
	   double threshold_scaling = fcal_config->FCAL_THRESHOLD_SCALING;
	   Ethreshold=FCAL_gain*integral_peak*MeV_FADC*(threshold*threshold_scaling - pedestal+gDRandom.SampleGaussian(pedestal_rms));
	   Erange *= FCAL_gain;
	   
	   if(fcal_config->FCAL_ADD_LIGHTGUIDE_HITS) {
	     hddm_s::FcalTruthLightGuideList lghits = titer->getFcalTruthLightGuides();
	     hddm_s::FcalTruthLightGuideList::iterator lgiter;
	     for (lgiter = lghits.begin(); lgiter != lghits.end(); lgiter++) {
	       E += lgiter->getE();
	     }
	   }
	   // Apply constant scale factor to MC energy. 06/22/2016 A. Subedi
	   E *= fcal_config->FCAL_MC_ESCALE;
	 }
	 else { // deal with insert
	   sigEfloor=fcal_config->INSERT_ENERGY_WIDTH_FLOOR;
	   sigEstat=fcal_config->INSERT_PHOT_STAT_COEF;
	 }
         
         double t = titer->getT(); 

         if(config->SMEAR_HITS) {
	   // Smear the timing and energy of the hit
	   t += gDRandom.SampleGaussian(fcal_config->FCAL_TSIGMA);
	   
	   // Energy width has stochastic and floor terms
	   double sigma = sqrt( pow(sigEstat,2)/titer->getE() 
			     + pow(sigEfloor,2));
	   E *= (1.0 + gDRandom.SampleGaussian(sigma));
	 }
	 	 
	 if (E > Erange)
	   E = Erange;
         
	 // Apply a single block threshold and energy range 
         // Scale threshold by gains
	 // Scale range by gains
	 if (E >= Ethreshold){
	   hddm_s::FcalHitList hits = iter->addFcalHits();
	   hits().setE(E);
	   hits().setT(t);
	 }
      }
      

      if (config->DROP_TRUTH_HITS)
         iter->deleteFcalTruthHits();
   }
   if (config->DROP_TRUTH_HITS) {
      hddm_s::ForwardEMcalList fcals = record->getForwardEMcals();
      if (fcals.size() > 0)
         fcals().deleteFcalTruthShowers();
   }
}
