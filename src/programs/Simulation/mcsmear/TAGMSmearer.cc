#include "TAGMSmearer.h"
#include "DANA/DEvent.h"


tagm_config_t *tagm_config_instance(0);

static int RANDOMIZE_ETAG_IN_EBIN(0);

//-----------
// tagm_config_t  (constructor)
//-----------
tagm_config_t::tagm_config_t(const std::shared_ptr<const JEvent>& event) {
	// default values
	TAGM_TSIGMA = 0.200;        // ns
	TAGM_FADC_TSIGMA = 0.350;   // ns
	TAGM_NPIX_PER_GEV = 1.e5;

    RANDOMIZE_ETAG_IN_EBIN = false;
    auto app = event->GetJApplication();
    app->SetDefaultParameter("TAGM:RANDOMIZE_ETAG_IN_EBIN",
                                RANDOMIZE_ETAG_IN_EBIN,
                                "Turn on/off randomization of tagged photon energy"
                                " in smeared tagger microscope hits."
                                " Set to \"1\" to turn on (it's off by default)");

   // enable on-the-fly bzip2 compression on output stream

    std::vector<std::map<std::string, double> > quality;
    if (DEvent::GetCalib(event, "/PHOTON_BEAM/microscope/fiber_quality", quality)) {
	   jout << "/PHOTON_BEAM/microscope/fiber_quality not used for this run" << endl;
    }
    for (int i=0; i < (int)quality.size(); ++i) {
       int rowcolumn = quality[i]["row"]*1000 + quality[i]["column"];
       fiber_quality[rowcolumn] = quality[i]["code"];
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
   std::vector<std::map<string, float> > micro_parms;
   DEvent::GetCalib(event, "PHOTON_BEAM/microscope/scaled_energy_range", micro_parms);
   for (unsigned int i=0; i < micro_parms.size(); ++i) {
      int col = micro_parms[i]["column"];
      double Emin = micro_parms[i]["xlow"] * endpoint_calib_GeV
                    + endpoint_energy_GeV - endpoint_calib_GeV;
      double Emax = micro_parms[i]["xhigh"] * endpoint_calib_GeV
                    + endpoint_energy_GeV - endpoint_calib_GeV;
      energy_range_GeV[col] = {Emin, Emax};
   }
   tagm_config_instance = this;
}

double TAGMSmearer::get_tagm_energy(int column, int mid_low_high_rand)
{
   // mid_low_high_rand: =0 returns the center energy of the bin,
   //                    =1 returns the low energy of the bin,
   //                    =2 returns the high energy of the bin,
   //                    =other returns a random value within the bin.
   if (tagm_config_instance == 0)
      return -1;
   else if (mid_low_high_rand == 0)
      return (tagm_config_instance->energy_range_GeV[column][0] +
              tagm_config_instance->energy_range_GeV[column][1]) / 2;
   else if (mid_low_high_rand == 1)
      return tagm_config_instance->energy_range_GeV[column][0];
   else if (mid_low_high_rand == 2)
      return tagm_config_instance->energy_range_GeV[column][1];
   else if (RANDOMIZE_ETAG_IN_EBIN != 0)
      return gDRandom.Uniform(tagm_config_instance->energy_range_GeV[column][0],
                              tagm_config_instance->energy_range_GeV[column][1]);
   return (tagm_config_instance->energy_range_GeV[column][0] +
           tagm_config_instance->energy_range_GeV[column][1]) / 2;
}

//-----------
// SmearEvent
//-----------
void TAGMSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::MicroChannelList tagms = record->getMicroChannels();
   hddm_s::MicroChannelList::iterator iter;
   for (iter = tagms.begin(); iter != tagms.end(); ++iter) {
      int column = iter->getColumn();
      int row = iter->getRow();
      iter->setE(get_tagm_energy(column));
      iter->deleteTaggerHits();
      hddm_s::TaggerTruthHitList thits = iter->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         if (tagm_config->fiber_quality[column + row*1000] != 1)
            continue;
         // smear the time
         double t = titer->getT();
         double tADC = titer->getT();
         double npe = titer->getDE() * tagm_config->TAGM_NPIX_PER_GEV;
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(tagm_config->TAGM_TSIGMA);
          	tADC += gDRandom.SampleGaussian(tagm_config->TAGM_FADC_TSIGMA);
          	npe = gDRandom.SamplePoisson(titer->getDE() * tagm_config->TAGM_NPIX_PER_GEV);
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
