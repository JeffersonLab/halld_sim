#include "FMWPCSmearer.h"

//-----------
// fmwpc_config_t  (constructor)
//-----------
fmwpc_config_t::fmwpc_config_t(const std::shared_ptr<const JEvent>& event) 
{
  // default values
  FMWPC_TSIGMA = 1.0;  // ns
  FMWPC_ASIGMA = 0.0; // ???
  FMWPC_THRESHOLD = 0.0; 
}



//-----------
// SmearEvent
//-----------
void FMWPCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::FmwpcChamberList chambers = record->getFmwpcChambers();
   hddm_s::FmwpcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {
      iter->deleteFmwpcHits();
      hddm_s::FmwpcTruthHitList thits = iter->getFmwpcTruthHits();
      hddm_s::FmwpcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
	hddm_s::FmwpcTruthHitQList &charges=titer->getFmwpcTruthHitQs();
         // smear the time and energy
         double t = titer->getT();
         double q = (charges.size()) ? charges.begin()->getQ() : 0.;
	 // Approximate drift time
	 const double v=0.0046; // drift velocity, cm/ns
	 double d = (charges.size()) ? charges.begin()->getD() : 0.;
	 double tdrift=d/v;
	 // Approximate longitudinal diffusion
	 double D=2.e-6; // cm^2/ns, based on 200 micron/cm for large drift time
	 double sigma_t=sqrt(2.*D*tdrift)/v;
	 t += tdrift+gDRandom.SampleGaussian(sigma_t);
         if(config->SMEAR_HITS) {
	   t += gDRandom.SampleGaussian(fmwpc_config->FMWPC_TSIGMA);
	   q += gDRandom.SampleGaussian(fmwpc_config->FMWPC_ASIGMA);
	 }
         if (q > fmwpc_config->FMWPC_THRESHOLD) {
	   hddm_s::FmwpcHitList hits = iter->addFmwpcHits();
	   hits().setT(t);
	   hits().setDE(0.); // Not used (SJT 2/28/22)
	   hddm_s::FmwpcHitQList charges=hits().addFmwpcHitQs(1);
	   charges(0).setQ(q);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteFmwpcTruthHits();
   }
}
