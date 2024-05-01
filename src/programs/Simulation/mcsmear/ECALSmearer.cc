
#include "ECALSmearer.h"
#include "DANA/DEvent.h"


DECALGeometry *ecalGeom = NULL;

//-----------
// ecal_config_t  (constructor)
//-----------
ecal_config_t::ecal_config_t(const std::shared_ptr<const JEvent>& event) {

  // Default Parameters

        ECAL_EN_SCALE  =  1.0962;
 	
	// Measured energy resolution
	ECAL_EN_P0     =  3.08e-2;
	ECAL_EN_P1     =  1.e-2;
	ECAL_EN_P2     =  0.7e-2;	

	// Energy deposition in Geant
	ECAL_EN_GP0     =  1.71216e-2;
	ECAL_EN_GP1     =  1.55070e-2;
	ECAL_EN_GP2     =  0.0;		
	
	// Time smearing factor
	ECAL_TSIGMA     =  0.4;
	
	
	// Single block energy threshold (applied after smearing)
	ECAL_BLOCK_THRESHOLD = 15.0*k_MeV;


        // Get values from CCDB

        cout << "Get ECAL/mc_energy parameters from CCDB..." << endl;

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

	cout<<"get ECAL/mc_time parameters from calibDB"<<endl;

	map<string, double> ecaltime;
	if(DEvent::GetCalib(event, "ECAL/mc_time", ecaltime)) {
	  jerr << "Problem loading ECAL/mc_time from CCDB!" << endl;
	} else {
	  ECAL_TSIGMA = ecaltime["ECAL_TSIGMA"];
	}
	
}



//-----------
// SmearEvemt
//-----------
void ECALSmearer::SmearEvent(hddm_s::HDDM *record){

#if 1 

  //   if (!ecalGeom)
  //   ecalGeom = new DECALGeometry();

  hddm_s::EcalBlockList blocks = record->getEcalBlocks();   
  hddm_s::EcalBlockList::iterator iter;
  for (iter = blocks.begin(); iter != blocks.end(); ++iter) {
    iter->deleteEcalHits();
    hddm_s::EcalTruthHitList thits = iter->getEcalTruthHits();   
    hddm_s::EcalTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      // Simulation simulates a grid of blocks for simplicity. 
      // Do not bother smearing inactive blocks. They will be
      // discarded in DEventSourceHDDM.cc while being read in
      // anyway.
      
      if (!ecalGeom->isBlockActive(iter->getRow(), iter->getColumn()))
		continue;

      
      // A.S.  new calibration of the ECAL
      double E = titer->getE();
      double t = titer->getT();
      
      E *= ecal_config->ECAL_EN_SCALE;
      
      if(config->SMEAR_HITS) {
	
	// Expected detector resolution
	double de_e_expect  =   pow(ecal_config->ECAL_EN_P0/sqrt(E),2) + 
	  pow(ecal_config->ECAL_EN_P1/E,2) + ecal_config->ECAL_EN_P2*ecal_config->ECAL_EN_P2; 
	
	// Subtract intrinsic Geant resolution
	double de_e_geant   =   pow(ecal_config->ECAL_EN_GP0/sqrt(E),2) + pow(ecal_config->ECAL_EN_GP1/E,2);
	
	double sig_res      =   sqrt(de_e_expect - de_e_geant);
	
	if(sig_res > 0) 
	  E *= (1. + gDRandom.SampleGaussian(sig_res));
	
	t += gDRandom.SampleGaussian(ecal_config->ECAL_TSIGMA);
	
      }
      
      
      
      // A.S.  Don't apply energy threshold at the moment

      //         if (E > ecal_config->ECAL_BLOCK_THRESHOLD) {
      hddm_s::EcalHitList hits = iter->addEcalHits();
      hits().setE(E);
      hits().setT(t);
      //         }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteEcalTruthHits();
  }
  if (config->DROP_TRUTH_HITS) {
    hddm_s::CrystalEcalList ecals = record->getCrystalEcals();   
    if (ecals.size() > 0)
      ecals().deleteEcalTruthShowers();
  }

#endif

}
