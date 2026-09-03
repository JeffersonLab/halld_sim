
#include "ECALSmearer.h"
#include "DANA/DEvent.h"


//-----------
// ecal_config_t  (constructor)
//-----------
ecal_config_t::ecal_config_t(const std::shared_ptr<const JEvent>& event, const DECALGeometry *ecalGeom) {

  // Default Parameters

        ECAL_EN_SCALE  =  1.0962;
 	
	// Measured energy resolution
	ECAL_EN_P0     =  3.08e-2;
	ECAL_EN_P1     =  1.e-2;
	ECAL_EN_P2     =  0.7e-2;	

	// Energy deposition in Geant
	ECAL_EN_GP0     =  1.71216e-2;
	ECAL_EN_GP1     =  1.e-2;
	ECAL_EN_GP2     =  0.0;		
	
	// Time smearing factor
	ECAL_TSIGMA     =  0.4;
	      
	// Single block energy threshold (applied after smearing)
	ECAL_ADC_THRESHOLD = 107;

	// Baseline fluctuation
	PED_SIGMA  = 0;

	ADC_EN_SCALE = 0.5447e-3;

       	INT_OVER_PEAK = 5.18;
	
        // Get values from CCDB

        cout << "Get ECAL/mc_energy parameters from CCDB ..." << endl;
	
	map<string, double> ecalparms;
	
	if(DEvent::GetCalib(event, "ECAL/mc_energy", ecalparms)) { 
	  jerr << "Problem loading ECAL/mc_energy from CCDB!" << endl;
	} else {
	  ECAL_EN_SCALE   = ecalparms["ECAL_EN_SCALE"]; 

	  ECAL_EN_P0    =  ecalparms["ECAL_EN_P0"]; 
	  ECAL_EN_P1    =  ecalparms["ECAL_EN_P1"]; 
	  ECAL_EN_P2    =  ecalparms["ECAL_EN_P2"]; 
	  
	  ECAL_EN_GP0   =  ecalparms["ECAL_EN_GP0"]; 
	  ECAL_EN_GP1   =  ecalparms["ECAL_EN_GP1"]; 
	  ECAL_EN_GP2   =  ecalparms["ECAL_EN_GP2"]; 
        }

	cout << "Get ECAL/mc_time parameters from CCDB ..." << endl;
	map<string, double> ecaltime;
	
	if(DEvent::GetCalib(event, "ECAL/mc_time", ecaltime)) {
	  jerr << "Problem loading ECAL/mc_time from CCDB!" << endl;
	} else {
	  ECAL_TSIGMA = ecaltime["ECAL_TSIGMA"];
	}
	
	cout << "Get ECAL/digi_scales parameters from CCDB ..." << endl;	
	map<string,double> scale_factors;
	
	if (DEvent::GetCalib(event, "/ECAL/digi_scales", scale_factors))
	  jout << "Error loading /ECAL/digi_scales !" << endl;
	if (scale_factors.find("ADC_EN_SCALE") != scale_factors.end())
	  ADC_EN_SCALE = scale_factors["ADC_EN_SCALE"];
	else
	  jerr << "Unable to get ADC_EN_SCALE from /ECAL/digi_scales !" << endl;
	
	cout << "Get ECAL/mc_parms from CCDB ..." << endl;	
	map<string, double> mcparms;

	if(DEvent::GetCalib(event, "ECAL/mc_parms", mcparms)) { 
	  jerr << "Problem loading ECAL/mc_parms from CCDB!" << endl;
	} else {
	  ECAL_ADC_THRESHOLD  = mcparms["THRESHOLD"]; 
	  INT_OVER_PEAK       = mcparms["INTEGRAL_PEAK"]; 
	  PED_SIGMA           = mcparms["PED_SIGMA"]; 
        }       
	
        int max_chan = DECALGeometry::kECALMaxChannels;

	cout << "Get ECAL/gains from CCDB ..." << endl;
	vector<double> ecal_gains_ch;
	
	if (DEvent::GetCalib(event, "/ECAL/gains", ecal_gains_ch)){
	  jout << "DECALHit_factory: Error loading /ECAL/gains !" << endl;
	  for (int ch = 0; ch < max_chan; ch ++) GAINS.push_back(1.);
	}
	else {
	  for (int ch = 0; ch < static_cast<int>(ecal_gains_ch.size()); ch++) {
	    GAINS.push_back(ecal_gains_ch[ch]);
	  }
	}
	
	cout << "Get ECAL/pedestals from CCDB ..." << endl;	
	vector<double> ecal_pedestals_ch;
	
	if (DEvent::GetCalib(event, "/ECAL/pedestals", ecal_pedestals_ch)){
	  jout << "DECALHit_factory: Error loading /ECAL/pedestals !" << endl;
	  for (int ch = 0; ch < max_chan; ch ++) PEDESTALS.push_back(100.);
	}
	else {
	  for (int ch = 0; ch < static_cast<int>(ecal_pedestals_ch.size()); ch++) {
	    PEDESTALS.push_back(ecal_pedestals_ch[ch]);
	  }
	}

	cout << "Get ECAL/bad_clock from CCDB ..." << endl;
	vector<double> ecal_bad_blocks_ch;
	
	if (DEvent::GetCalib(event, "/ECAL/bad_block", ecal_bad_blocks_ch)){
	  jout << "DECALHit_factory: Error loading /ECAL/bad_block !" << endl;
	  for (int ch = 0; ch < max_chan; ch ++) BAD_BLOCKS.push_back(0.);
	}
	else {
	  for (int ch = 0; ch < static_cast<int>(ecal_bad_blocks_ch.size()); ch++) {
	    BAD_BLOCKS.push_back(ecal_bad_blocks_ch[ch]);
	  }
	}	
	
}



//-----------
// SmearEvent
//-----------
void ECALSmearer::SmearEvent(hddm_s::HDDM *record){

  hddm_s::EcalBlockList blocks = record->getEcalBlocks();   
  hddm_s::EcalBlockList::iterator iter;
  for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
    iter->deleteEcalHits();
    hddm_s::EcalTruthHitList thits = iter->getEcalTruthHits();   
    hddm_s::EcalTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      
      // A.S.  new calibration of the ECAL
      double E = titer->getE();
      double t = titer->getT();

      int column=iter->getColumn();
      int row=iter->getRow();
      
      int chan = ecalGeom->channel(row, column);

      double en_scale_cor = ecal_config->ECAL_EN_SCALE;
      
      E *= en_scale_cor;
      
      if(config->SMEAR_HITS) {
	
	// Expected detector resolution
	double de_e_expect  =   pow(ecal_config->ECAL_EN_P0/sqrt(E),2) + 
	  pow(ecal_config->ECAL_EN_P1/E,2) + ecal_config->ECAL_EN_P2*ecal_config->ECAL_EN_P2; 
	
	// Subtract intrinsic Geant resolution
	double de_e_geant   =   pow(ecal_config->ECAL_EN_GP0/sqrt(E),2) + pow(ecal_config->ECAL_EN_GP1/E,2);
	
	double sig_res = 0;
	
	if((de_e_expect - de_e_geant) > 0)
	  sig_res      =   sqrt(de_e_expect - de_e_geant);
	
	if(sig_res > 0) 
	  E *= (1. + gDRandom.SampleGaussian(sig_res));
	
	t += gDRandom.SampleGaussian(ecal_config->ECAL_TSIGMA);	
      }
      
      
      // A.S.  Calculate energy threshold
      int bad_block           =   ecal_config->BAD_BLOCKS.at(chan);
      double adc_threshold    =   ecal_config->ECAL_ADC_THRESHOLD;
      double pedestal         =   ecal_config->PEDESTALS.at(chan);
      double pedestal_sigma   =   ecal_config->PED_SIGMA;
      double gain             =   ecal_config->GAINS.at(chan);
      double adc_en_scale     =   ecal_config->ADC_EN_SCALE;
      double intOverPeak      =   ecal_config->INT_OVER_PEAK;
      double pedestal_fluct   =   gDRandom.SampleGaussian(pedestal_sigma);
      
      double threshold = adc_en_scale*intOverPeak*en_scale_cor*gain*(adc_threshold - pedestal + pedestal_fluct);
      
      double baseline_shift = adc_threshold - pedestal;
      
      if(bad_block > 0) continue;       // Bad block

      // Check the energy threshold. Do not produce hits if the baseline is above the threshold. 
      if ((E > threshold) && (baseline_shift > 0)) {    // Check block threshold and baseline
	hddm_s::EcalHitList hits = iter->addEcalHits();
	hits().setE(E);
	hits().setT(t);
	
      }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteEcalTruthHits();
  }
  if (config->DROP_TRUTH_HITS) {
    hddm_s::CrystalEcalList ecals = record->getCrystalEcals();   
    if (ecals.size() > 0)
      ecals().deleteEcalTruthShowers();
  }
  
}
