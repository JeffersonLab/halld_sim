#include "DIRCSmearer.h"

//-----------
// dirc_config_t  (constructor)
//-----------
dirc_config_t::dirc_config_t(JEventLoop *loop, DDIRCGeometry *dircGeom)
{
        // default values
        DIRC_TSIGMA           = 0.5; // 0.5 ns 

	// Get values from CCDB
	cout<<"get DIRC/mc_timing_smear parameters from calibDB"<<endl;
	map<string, double> dircmctimingsmear;
	if(loop->GetCalib("DIRC/mc_timing_smear", dircmctimingsmear)) {
		jerr << "Problem loading DIRC/mc_timing_smear from CCDB!" << endl;
	} else {
		DIRC_TSIGMA = dircmctimingsmear["DIRC_TSIGMA"];
	}

	// get per-pixel efficiencies from CCDB
}

//-----------
// SmearEvent
//-----------
void DIRCSmearer::SmearEvent(hddm_s::HDDM *record)
{
	hddm_s::DircTruthPmtHitList truthPmtHits = record->getDircTruthPmtHits();
	hddm_s::DircTruthPmtHitList::iterator iter;
	for (iter = truthPmtHits.begin(); iter != truthPmtHits.end(); ++iter) {
		iter->deleteDircPmtHits();
		
		// add per-pixel efficiencies from MAPMT test data
		//if (config->APPLY_EFFICIENCY_CORRECTIONS && !gDRandom.DecideToAcceptHit(dirc_config->GetEfficiencyCorrectionFactor(iter->getCh())) ) {
		//	continue;
		//}

		double t = iter->getT();
		double t_fixed = iter->getT_fixed();
		int ch = iter->getCh();
		
		if(config->SMEAR_HITS) {
			// Smear the timing of the hit
                        t += gDRandom.SampleGaussian(dirc_config->DIRC_TSIGMA);
			
			// Add cross talk here?
		}
		
		hddm_s::DircPmtHitList hits = iter->addDircPmtHits();
		hits().setT(t);
		hits().setCh(ch);
	}
}
