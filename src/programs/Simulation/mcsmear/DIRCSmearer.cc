#include "DIRCSmearer.h"
#include "DANA/DEvent.h"

//-----------
// dirc_config_t  (constructor)
//-----------
dirc_config_t::dirc_config_t(const std::shared_ptr<const JEvent>& event)
{
        // default values
        DIRC_TSIGMA           = 0.5; // 0.5 ns 
	DIRC_MAX_CHANNELS     = 108*64; 

#if 0
	// Get values from CCDB
	cout<<"get DIRC/mc_timing_smear parameters from calibDB"<<endl;
	map<string, double> dircmctimingsmear;
	if(DEvent::GetCalib(event, "DIRC/mc_timing_smear", dircmctimingsmear)) {
		jerr << "Problem loading DIRC/mc_timing_smear from CCDB!" << endl;
	} else {
		DIRC_TSIGMA = dircmctimingsmear["DIRC_TSIGMA"];
	}
#endif

	// get DIRC channel status from DB
	vector<int> new_status(DIRC_MAX_CHANNELS);
	dChannelStatus.push_back(new_status); 
	dChannelStatus.push_back(new_status);
	if (DEvent::GetCalib(event, "/DIRC/North/channel_status", dChannelStatus[0]))
		jout << "Error loading /DIRC/North/channel_status !" << endl;
	if (DEvent::GetCalib(event, "/DIRC/South/channel_status", dChannelStatus[1]))
		jout << "Error loading /DIRC/South/channel_status !" << endl;
	
	// get per-pixel efficiencies from CCDB
}

//-----------
// SmearEvent
//-----------
void DIRCSmearer::SmearEvent(hddm_s::HDDM *record)
{
#ifdef SMEARDIRC
	hddm_s::DIRCList dirc = record->getDIRCs();
	if(dirc.size() > 0) 
		dirc().deleteDircPmtHits();

	hddm_s::DircTruthPmtHitList truthPmtHits = record->getDircTruthPmtHits();
	hddm_s::DircTruthPmtHitList::iterator iter;
	for (iter = truthPmtHits.begin(); iter != truthPmtHits.end(); ++iter) {
		
		// add per-pixel efficiencies from MAPMT test data
		//if (config->APPLY_EFFICIENCY_CORRECTIONS && !gDRandom.DecideToAcceptHit(dirc_config->GetEfficiencyCorrectionFactor(iter->getCh())) ) {
		//	continue;
		//}
		
		double t = iter->getT();
		int ch = iter->getCh();
		
		if(config->SMEAR_HITS) {
			// Smear the timing of the hit
			t += gDRandom.SampleGaussian(dirc_config->DIRC_TSIGMA);
			
			// Add cross talk here?
			
			// Remove pixels with bad status
			int box = (iter->getCh() < dirc_config->DIRC_MAX_CHANNELS) ? 1 : 0;
			int channel = iter->getCh() % dirc_config->DIRC_MAX_CHANNELS;
			dirc_status_state status = static_cast<dirc_status_state>(dirc_config->dChannelStatus[box][channel]);
			if ( (status==BAD) || (status==NOISY) ) {
				continue;
			}
		}
		
		hddm_s::DircPmtHitList hits = dirc().addDircPmtHits();
		hits().setT(t);
		hits().setCh(ch);
	}
#endif
    if (config->DROP_TRUTH_HITS) {
	    hddm_s::DIRCList dircs = record->getDIRCs();
        if  (dircs.size() > 0) {
            dircs().deleteDircTruthPmtHits();
            dircs().deleteDircTruthBarHits();
        }
    }
}
