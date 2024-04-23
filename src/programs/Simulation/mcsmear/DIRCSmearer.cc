#include "DIRCSmearer.h"

//-----------
// dirc_config_t  (constructor)
//-----------
dirc_config_t::dirc_config_t(JEventLoop *loop)
{
        // default values
	DIRC_EFFIC_SCALE      = 1.5;
        DIRC_TSIGMA           = 0.7; // 0.7 ns
	DIRC_EFFIC_DEGRADE    = 1.0;
	DIRC_MAX_CHANNELS     = 108*64; 

	// get time smearing from CCDB
	cout<<"get DIRC/mc_parms parameters from calibDB"<<endl;
	map<string, double> mc_parms;
	if(loop->GetCalib("DIRC/mc_parms", mc_parms)) {
		jerr << "Problem loading DIRC/mc_parms from CCDB!" << endl;
	} else {
		DIRC_EFFIC_SCALE = mc_parms["PAR0"];
		DIRC_TSIGMA = mc_parms["PAR1"];
		//DIRC_EFFIC_DEGRADE = mc_parms["PAR2"];
	}

	// get DIRC channel status and efficiency from DB
	vector<int> new_status(DIRC_MAX_CHANNELS);
	vector<float> new_effic(DIRC_MAX_CHANNELS);
	dChannelStatus.push_back(new_status); 
	dChannelStatus.push_back(new_status);
	dChannelEffic.push_back(new_effic);
        dChannelEffic.push_back(new_effic);
	if (loop->GetCalib("/DIRC/North/channel_status", dChannelStatus[0]))
		jout << "Error loading /DIRC/North/channel_status !" << endl;
	if (loop->GetCalib("/DIRC/South/channel_status", dChannelStatus[1]))
		jout << "Error loading /DIRC/South/channel_status !" << endl;
	if (loop->GetCalib("/DIRC/North/channel_effic", dChannelEffic[0]))
                jout << "Error loading /DIRC/North/channel_effic !" << endl;
        if (loop->GetCalib("/DIRC/South/channel_effic", dChannelEffic[1]))
                jout << "Error loading /DIRC/South/channel_effic !" << endl;
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
		
		double t = iter->getT();
		int ch = iter->getCh();
		
		if(config->SMEAR_HITS) {
			// Smear the timing of the hit
			t += gDRandom.SampleGaussian(dirc_config->DIRC_TSIGMA);

			// Remove pixels with bad status
			int box = (iter->getCh() < dirc_config->DIRC_MAX_CHANNELS) ? 1 : 0;
			int channel = iter->getCh() % dirc_config->DIRC_MAX_CHANNELS;
			dirc_status_state status = static_cast<dirc_status_state>(dirc_config->dChannelStatus[box][channel]);
			if ( (status==BAD) || (status==NOISY) ) {
				continue;
			}

			// Add per-pixel efficiencies from MAPMT test data
			if (config->APPLY_EFFICIENCY_CORRECTIONS && !gDRandom.DecideToAcceptHit(dirc_config->dChannelEffic[box][channel] / dirc_config->DIRC_EFFIC_SCALE * dirc_config->DIRC_EFFIC_DEGRADE)) {
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
