
#include "CCALSmearer.h"

DCCALGeometry *ccalGeom = NULL;

//-----------
// ccal_config_t  (constructor)
//-----------
ccal_config_t::ccal_config_t(JEventLoop *loop) {
	
	// Default Parameters
	
	CCAL_EN_SCALE = 1.0962;
	
	// Measured energy resolution
	CCAL_EN_P0    = 3.08e-2;
	CCAL_EN_P1    = 1.e-2;
	CCAL_EN_P2    = 0.7e-2;
	
	// Energy deposition in Geant
	CCAL_EN_GP0   = 1.71216e-2;
	CCAL_EN_GP1   = 1.55070e-2;
	CCAL_EN_GP2   = 0.0;
	
	// Time smearing factor
	CCAL_TSIGMA   = 0.2;
	
	// Single block energy threshold (applied after smearing)
	CCAL_BLOCK_THRESHOLD = 15.0*k_MeV;
	
	// Get values from CCDB
	
	cout << "Get CCAL/mc_energy parameters from CCDB..." << endl;
	
	map<string, double> ccalparms;
	
	if(loop->GetCalib("CCAL/mc_energy", ccalparms)) {
		jerr << "Problem loading CCAL/mc_energy from CCDB!" << endl;
	} else {
		CCAL_EN_SCALE = ccalparms["CCAL_EN_SCALE"];
		
		CCAL_EN_P0    = ccalparms["CCAL_EN_P0"];
		CCAL_EN_P1    = ccalparms["CCAL_EN_P1"];
		CCAL_EN_P2    = ccalparms["CCAL_EN_P2"];
		
		CCAL_EN_GP0   = ccalparms["CCAL_EN_GP0"];
		CCAL_EN_GP1   = ccalparms["CCAL_EN_GP1"];
		CCAL_EN_GP2   = ccalparms["CCAL_EN_GP2"];
	}
	
	cout<<"get CCAL/mc_time parameters from calibDB"<<endl;
	
	map<string, double> ccaltime;
	if(loop->GetCalib("CCAL/mc_time", ccaltime)) {
		jerr << "Problem loading CCAL/mc_time from CCDB!" << endl;
	} else {
		CCAL_TSIGMA = ccaltime["CCAL_TSIGMA"];
	}
	
	cout<<"get CCAL/gains from calibDB"<<endl;
	vector <double> CCAL_GAINS_TEMP;
	if(loop->GetCalib("CCAL/gains", CCAL_GAINS_TEMP)) {
		jerr << "Problem loading CCAL/gains from CCDB!" << endl;
	} else {
		for(unsigned int i = 0; i < CCAL_GAINS_TEMP.size(); i++) {
			CCAL_GAINS.push_back(CCAL_GAINS_TEMP.at(i));
		}
	}
	
	cout<<"get CCAL/pedestals from calibDB"<<endl;
	vector <double> CCAL_PEDS_TEMP;
	if(loop->GetCalib("CCAL/pedestals", CCAL_PEDS_TEMP)) {
		jerr << "Problem loading CCAL/pedestals from CCDB!" << endl;
	} else {
		for(unsigned int i = 0; i < CCAL_PEDS_TEMP.size(); i++) {
			CCAL_PEDS.push_back(CCAL_PEDS_TEMP.at(i));
		}
	}
	
	// Hard code these values for now. Eventually they should be stored and accessed from the database:
	CCAL_INTEGRAL_PEAK = 5.7;
	CCAL_THRESHOLD     = 108.;
	CCAL_PED_RMS       = 1.0;
	
	cout<<"get CCAL/digi_scales from calibDB"<<endl;
	map<string, double> ccaldigiscales;
	if(loop->GetCalib("CCAL/digi_scales", ccaldigiscales)) {
		jerr << "Problem loading CCAL/digi_scales from CCDB!" << endl;
	} else {
		for(unsigned int i = 0; i < CCAL_PEDS_TEMP.size(); i++) {
			CCAL_ADC_ASCALE = ccaldigiscales["ADC_EN_SCALE"];
		}
	}
}


//-----------
// SmearEvent
//-----------
void CCALSmearer::SmearEvent(hddm_s::HDDM *record){
	
	//if (!ccalGeom)
	//ccalGeom = new DCCALGeometry();
	
	hddm_s::CcalBlockList blocks = record->getCcalBlocks();
	hddm_s::CcalBlockList::iterator iter;
	for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
		iter->deleteCcalHits();
		hddm_s::CcalTruthHitList thits = iter->getCcalTruthHits();
		hddm_s::CcalTruthHitList::iterator titer;
		
		int row    = iter->getRow();
		int column = iter->getColumn();
		
		for (titer = thits.begin(); titer != thits.end(); ++titer) {
			// Simulation simulates a grid of blocks for simplicity.
			// Do not bother smearing inactive blocks. They will be
			// discarded in DEventSourceHDDM.cc while being read in
			// anyway.
			
			if (!ccalGeom->isBlockActive(iter->getRow(), iter->getColumn()))
				continue;
			
			// A.S.  new calibration of the CCAL
			double E = titer->getE();
			double t = titer->getT();
			
			E *= ccal_config->CCAL_EN_SCALE;
			
			if(config->SMEAR_HITS) {
				
				// Expected detector resolution
				double de_e_expect  =   pow(ccal_config->CCAL_EN_P0/sqrt(E),2) + 
					pow(ccal_config->CCAL_EN_P1/E,2) + ccal_config->CCAL_EN_P2*ccal_config->CCAL_EN_P2;
				
				// Subtract intrinsic Geant resolution
				double de_e_geant   =   pow(ccal_config->CCAL_EN_GP0/sqrt(E),2) + pow(ccal_config->CCAL_EN_GP1/E,2);
				
				double sig_res      =   sqrt(de_e_expect - de_e_geant);
				
				if(sig_res > 0)
					E *= (1. + gDRandom.SampleGaussian(sig_res));
				
				t += gDRandom.SampleGaussian(ccal_config->CCAL_TSIGMA);
			}
			
			// D.Smith (7/18/2024): Apply channel-dependent energy thresholds:
			
			double integral_peak_ratio = ccal_config->CCAL_INTEGRAL_PEAK;
			double threshold           = ccal_config->CCAL_THRESHOLD;
			double MeV_FADC            = ccal_config->CCAL_ADC_ASCALE;
			double pedestal_rms        = ccal_config->CCAL_PED_RMS;
			
			int channelnum = ccalGeom->channel(row, column);
			double loc_gain = ccal_config->CCAL_GAINS[channelnum];
			double loc_ped  = ccal_config->CCAL_PEDS[channelnum];
			
			double Ethreshold = loc_gain * integral_peak_ratio * MeV_FADC * (threshold - loc_ped 
				+ gDRandom.SampleGaussian(pedestal_rms)); // in units of MeV
			
			// A.S.  Don't apply energy threshold at the moment
			
			if(E >= (Ethreshold*1.e-3)) {
				hddm_s::CcalHitList hits = iter->addCcalHits();
				hits().setE(E*1000.);
				hits().setT(t);
			}
		}
		
		if (config->DROP_TRUTH_HITS)
			iter->deleteCcalTruthHits();
	}
	if (config->DROP_TRUTH_HITS) {
		hddm_s::ComptonEMcalList ccals = record->getComptonEMcals();
		if (ccals.size() > 0)
			ccals().deleteCcalTruthShowers();
	}
}
