
#include "CCALSmearer.h"
#include "DANA/DEvent.h"


DCCALGeometry *ccalGeom = NULL;

//-----------
// ccal_config_t  (constructor)
//-----------
ccal_config_t::ccal_config_t(const std::shared_ptr<const JEvent>& event) {

  // Default Parameters

        CCAL_EN_SCALE  =  1.0962;
 	
	// Measured energy resolution
	CCAL_EN_P0     =  3.08e-2;
	CCAL_EN_P1     =  1.e-2;
	CCAL_EN_P2     =  0.7e-2;	

	// Energy deposition in Geant
	CCAL_EN_GP0     =  1.71216e-2;
	CCAL_EN_GP1     =  1.55070e-2;
	CCAL_EN_GP2     =  0.0;		
	
	// Time smearing factor
	CCAL_TSIGMA     =  0.4;
	
	
	// Single block energy threshold (applied after smearing)
	CCAL_BLOCK_THRESHOLD = 15.0*k_MeV;



        // Get values from CCDB

        cout << "Get CCAL/mc_energy parameters from CCDB..." << endl;

	map<string, double> ccalparms;

	if(DEvent::GetCalib(event, "CCAL/mc_energy", ccalparms)) { 
	  jerr << "Problem loading CCAL/mc_energy from CCDB!" << endl;
	} else {
	  CCAL_EN_SCALE   = ccalparms["CCAL_EN_SCALE"]; 

	  CCAL_EN_P0    =  ccalparms["CCAL_EN_P0"]; 
	  CCAL_EN_P1    =  ccalparms["CCAL_EN_P1"]; 
	  CCAL_EN_P2    =  ccalparms["CCAL_EN_P2"]; 

	  CCAL_EN_GP0   =  ccalparms["CCAL_EN_GP0"]; 
	  CCAL_EN_GP1   =  ccalparms["CCAL_EN_GP1"]; 
	  CCAL_EN_GP2   =  ccalparms["CCAL_EN_GP2"]; 
        }

	cout<<"get CCAL/mc_time parameters from calibDB"<<endl;

	map<string, double> ccaltime;
	if(DEvent::GetCalib(event, "CCAL/mc_time", ccaltime)) {
	  jerr << "Problem loading CCAL/mc_time from CCDB!" << endl;
	} else {
	  CCAL_TSIGMA = ccaltime["CCAL_TSIGMA"];
	}
	

}



//-----------
// SmearEvemt
//-----------
void CCALSmearer::SmearEvent(hddm_s::HDDM *record){
  
  //   if (!ccalGeom)
  //   ccalGeom = new DCCALGeometry();
  

  hddm_s::CcalBlockList blocks = record->getCcalBlocks();   
  hddm_s::CcalBlockList::iterator iter;
  for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
    iter->deleteCcalHits();
    hddm_s::CcalTruthHitList thits = iter->getCcalTruthHits();   
    hddm_s::CcalTruthHitList::iterator titer;
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
      
      
      
      // A.S.  Don't apply energy threshold at the moment

      //         if (E > ccal_config->CCAL_BLOCK_THRESHOLD) {
      hddm_s::CcalHitList hits = iter->addCcalHits();
      hits().setE(E*1000.);
      hits().setT(t);
      //         }
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
