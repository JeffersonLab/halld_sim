#include "TAGHSmearer.h"
#include "DANA/DEvent.h"

tagh_config_t *tagh_config_instance(0);

static int RANDOMIZE_ETAG_IN_EBIN(0);

//-----------
// tagh_config_t  (constructor)
//-----------
tagh_config_t::tagh_config_t(const std::shared_ptr<const JEvent>& event) 
{
	// default values
	TAGH_TSIGMA = 0.350;        // ns
	TAGH_FADC_TSIGMA = 0.450;   // ns
	TAGH_NPE_PER_GEV = 5.e5;

    RANDOMIZE_ETAG_IN_EBIN = false;
    auto app = event->GetJApplication();
    app->SetDefaultParameter("TAGH:RANDOMIZE_ETAG_IN_EBIN",
                                RANDOMIZE_ETAG_IN_EBIN,
                                "Turn on/off randomization of tagged photon energy"
                                " within the given energy bin for hits in the hodoscope."
                                " Set to \"1\" to turn on (it's off by default)");

   // enable on-the-fly bzip2 compression on output stream
    std::vector<std::map<std::string, double> > quality;
    if (DEvent::GetCalib(event, "/PHOTON_BEAM/hodoscope/counter_quality", quality)) {
	   jout << "/PHOTON_BEAM/hodoscope/counter_quality not used for this run" << endl;
    }
    for (int i=0; i < (int)quality.size(); ++i) {
       int id = quality[i]["id"];
       counter_quality[id] = quality[i]["code"];
    }

   std::map<string, float> beam_parms;
   DEvent::GetCalib(event, "PHOTON_BEAM/endpoint_energy", beam_parms);
   endpoint_energy_GeV = beam_parms.at("PHOTON_BEAM_ENDPOINT_ENERGY");
   std::map<string, float> beam_calib;
   DEvent::GetCalib(event, "PHOTON_BEAM/hodoscope/endpoint_calib", beam_calib);
   endpoint_calib_GeV = endpoint_energy_GeV;
   if (beam_calib.find("TAGGER_CALIB_ENERGY") != beam_calib.end()) {
      endpoint_calib_GeV = beam_calib.at("TAGGER_CALIB_ENERGY");
   }
   std::vector<std::map<string, float> > hodo_parms;
   DEvent::GetCalib(event, "PHOTON_BEAM/hodoscope/scaled_energy_range", hodo_parms);
   for (unsigned int i=0; i < hodo_parms.size(); ++i) {
      int col = hodo_parms[i]["counter"];
      double Emin = hodo_parms[i]["xlow"] * endpoint_calib_GeV
                    + endpoint_energy_GeV - endpoint_calib_GeV;
      double Emax = hodo_parms[i]["xhigh"] * endpoint_calib_GeV
                    + endpoint_energy_GeV - endpoint_calib_GeV;
      energy_range_GeV[col] = {Emin, Emax};
   }
   tagh_config_instance = this;
}

double TAGHSmearer::get_tagh_energy(int counter, int mid_low_high_rand)
{
   // mid_low_high_rand: =0 returns the center energy of the bin,
   //                    =1 returns the low energy of the bin,
   //                    =2 returns the high energy of the bin,
   //                    =other returns a random value within the bin.
   if (tagh_config_instance == 0)
      return -1;
   else if (mid_low_high_rand == 0)
      return (tagh_config_instance->energy_range_GeV[counter][0] +
              tagh_config_instance->energy_range_GeV[counter][1]) / 2;
   else if (mid_low_high_rand == 1)
      return tagh_config_instance->energy_range_GeV[counter][0];
   else if (mid_low_high_rand == 2)
      return tagh_config_instance->energy_range_GeV[counter][1];
   else if (RANDOMIZE_ETAG_IN_EBIN != 0)
      return gDRandom.Uniform(tagh_config_instance->energy_range_GeV[counter][0],
                              tagh_config_instance->energy_range_GeV[counter][1]);
   return (tagh_config_instance->energy_range_GeV[counter][0] +
           tagh_config_instance->energy_range_GeV[counter][1]) / 2;
}

//-----------
// SmearEvent
//-----------
void TAGHSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::HodoChannelList taghs = record->getHodoChannels();
   hddm_s::HodoChannelList::iterator iter;
   for (iter = taghs.begin(); iter != taghs.end(); ++iter) {
      int counter = iter->getCounterId();
      iter->setE(get_tagh_energy(counter));
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
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
