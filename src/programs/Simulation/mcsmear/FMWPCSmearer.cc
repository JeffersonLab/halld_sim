#include "FMWPCSmearer.h"

//-----------
// fmwpc_config_t  (constructor)
//-----------
fmwpc_config_t::fmwpc_config_t(JEventLoop *loop) 
{
  // default values
  FMWPC_TSIGMA = 1.0;  // ns
  FMWPC_ASIGMA = 0.0; // ???
  FMWPC_THRESHOLD = 0.0;
  FMWPC_INTEGRAL_TO_AMPLITUDE=1./28.8; // copied from CDC
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
         // smear the time and energy
         double t = titer->getT();
         double q = titer->getQ();
	 // Approximate drift time
	 const double v=0.0046; // drift velocity, cm/ns
	 double tdrift=titer->getD()/v;
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
	   hits().setQ(q);
	   double amp=q*fmwpc_config->FMWPC_INTEGRAL_TO_AMPLITUDE;
	   hits().setAmp(amp);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteFmwpcTruthHits();
   }
}
