#include "GEMTRDSmearer.h"

//-----------
// gemtrd_config_t  (constructor)
//-----------
gemtrd_config_t::gemtrd_config_t(JEventLoop *loop) 
{
  // default values
  GEMTRD_TSIGMA = 1.0;  // ns
  GEMTRD_ASIGMA = 0.0; // ???
  GEMTRD_XYSIGMA = 0.01; // cm
  GEMTRD_THRESHOLD = 0.0;
  GEMTRD_INTEGRAL_TO_AMPLITUDE=1./28.8; // copied from CDC
}



//-----------
// SmearEvent
//-----------
void GEMTRDSmearer::SmearEvent(hddm_s::HDDM *record)
{
  const double v=0.0033; // drift velocity, cm/ns
  // from figure 8a in arXiv:1110.6761 at E=1.5 keV/cm
  
  const double D=5.2e-7; // longitudinal diffusion coefficient, cm^2/ns
  // based on 125 micron/cm longitinudal diffusion from figure 12a 
  // in arXiv:1110.6761 at E=1.5 keV/cm
  
  const double Dt=4.6e-6; // ?? transverse diffusion coefficient, cm^2/ns
  // based on 375 micron/cm transverse diffusion from figure 16a
  // in arXiv:1110.6761 at E=1.5 keV/cm
  
  hddm_s::GemtrdChamberList chambers = record->getGemtrdChambers();
  hddm_s::GemtrdChamberList::iterator iter;
  for (iter = chambers.begin(); iter != chambers.end(); ++iter) {
    iter->deleteGemtrdHits();
    hddm_s::GemtrdTruthHitList thits = iter->getGemtrdTruthHits();
    hddm_s::GemtrdTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      // smear the time and energy
      double t = titer->getT();
      double q = titer->getQ();
      double x = titer->getX();
      double y = titer->getY();
      
      // Approximate drift time
      double tdrift=titer->getD()/v;
      t += tdrift;
      
      if(config->SMEAR_HITS) {
	t += gDRandom.SampleGaussian(gemtrd_config->GEMTRD_TSIGMA);

	// Approximate longitudinal diffusion
	double sigma_t=(tdrift>0.) ? sqrt(2.*D*tdrift)/v : 0.;
	t+=gDRandom.SampleGaussian(sigma_t);
	
	x += gDRandom.SampleGaussian(gemtrd_config->GEMTRD_XYSIGMA);
	y += gDRandom.SampleGaussian(gemtrd_config->GEMTRD_XYSIGMA);
	
	// Approximate transverse diffusion
	double sigma_xy=(tdrift>0) ? sqrt(2.*Dt*tdrift) : 0.;
	x += gDRandom.SampleGaussian(sigma_xy);
	y += gDRandom.SampleGaussian(sigma_xy);

	q += gDRandom.SampleGaussian(gemtrd_config->GEMTRD_ASIGMA);
      }
      if (q > gemtrd_config->GEMTRD_THRESHOLD) {
	hddm_s::GemtrdHitList hits = iter->addGemtrdHits();
	hits().setT(t);
	hits().setQ(q);
	hits().setX(x);
	hits().setY(y);
      }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteGemtrdTruthHits();
  }
}
