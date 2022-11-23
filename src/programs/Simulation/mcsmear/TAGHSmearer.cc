#include "TAGHSmearer.h"


//-----------
// tagh_config_t  (constructor)
//-----------
tagh_config_t::tagh_config_t(JEventLoop *loop) 
{
	// default values
	TAGH_TSIGMA = 0.350;        // ns
	TAGH_FADC_TSIGMA = 0.450;   // ns
	TAGH_NPE_PER_GEV = 5.e5;

    std::vector<std::map<std::string, double> > quality;
    if (loop->GetCalib("/PHOTON_BEAM/hodoscope/counter_quality", quality)) {
	   jout << "/PHOTON_BEAM/hodoscope/counter_quality not used for this run" << endl;
    }
    for (int i=0; i < (int)quality.size(); ++i) {
       int id = quality[i]["id"];
       counter_quality[id] = quality[i]["code"];
    }
}


//-----------
// SmearEvent
//-----------
void TAGHSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::HodoChannelList taghs = record->getHodoChannels();
   hddm_s::HodoChannelList::iterator iter;
   for (iter = taghs.begin(); iter != taghs.end(); ++iter) {
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         int counter = *(int*)titer->getAttribute("counterId");
         if (tagh_config->counter_quality[counter] != 1)
            continue;
         // smear the time
         double t = titer->getT();
         double tADC = titer->getT();
         double npe = titer->getDE() * tagh_config->TAGH_NPE_PER_GEV;

         if(config->SMEAR_HITS) {
        	t += gDRandom.SampleGaussian(tagh_config->TAGH_TSIGMA);
         	tADC += gDRandom.SampleGaussian(tagh_config->TAGH_FADC_TSIGMA);
         	npe = gDRandom.SamplePoisson(titer->getDE() * tagh_config->TAGH_NPE_PER_GEV);
		 }
         hddm_s::TaggerHitList hits = iter->addTaggerHits();
         hits().setT(t);
         hits().setTADC(tADC);
         hits().setNpe(npe);
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteTaggerTruthHits();
   }
}
