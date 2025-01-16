#include "TRDSmearer.h"

//-----------
// trd_config_t  (constructor)
//-----------
trd_config_t::trd_config_t(JEventLoop *loop) 
{
  // default values
  TRD_TSIGMA = 1.0;  // ns
  TRD_ASIGMA = 0.0; // ???
  TRD_XYSIGMA = 0.01; // cm
  TRD_THRESHOLD = 0.0;
  TRD_INTEGRAL_TO_AMPLITUDE=1./28.8; // copied from CDC
}



//-----------
// SmearEvent
//-----------
void TRDSmearer::SmearEvent(hddm_s::HDDM *record)
{
  const double v=0.0033; // drift velocity, cm/ns
  // from figure 8a in arXiv:1110.6761 at E=1.5 keV/cm
  
  const double D=5.2e-7; // longitudinal diffusion coefficient, cm^2/ns
  // based on 125 micron/cm longitinudal diffusion from figure 12a 
  // in arXiv:1110.6761 at E=1.5 keV/cm
  
  const double Dt=4.6e-6; // ?? transverse diffusion coefficient, cm^2/ns
  // based on 375 micron/cm transverse diffusion from figure 16a
  // in arXiv:1110.6761 at E=1.5 keV/cm
  
  hddm_s::TrdChamberList chambers = record->getTrdChambers();
  hddm_s::TrdChamberList::iterator iter;
  for (iter = chambers.begin(); iter != chambers.end(); ++iter) {
    iter->deleteTrdHits();
    hddm_s::TrdTruthHitList thits = iter->getTrdTruthHits();
    hddm_s::TrdTruthHitList::iterator titer;
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
	t += gDRandom.SampleGaussian(trd_config->TRD_TSIGMA);

	// Approximate longitudinal diffusion
	double sigma_t=(tdrift>0.) ? sqrt(2.*D*tdrift)/v : 0.;
	t+=gDRandom.SampleGaussian(sigma_t);
	
	x += gDRandom.SampleGaussian(trd_config->TRD_XYSIGMA);
	y += gDRandom.SampleGaussian(trd_config->TRD_XYSIGMA);
	
	// Approximate transverse diffusion
	double sigma_xy=(tdrift>0) ? sqrt(2.*Dt*tdrift) : 0.;
	x += gDRandom.SampleGaussian(sigma_xy);
	y += gDRandom.SampleGaussian(sigma_xy);

	q += gDRandom.SampleGaussian(trd_config->TRD_ASIGMA);
      }
      if (q > trd_config->TRD_THRESHOLD) {
	hddm_s::TrdHitList xhits = iter->addTrdHits();
	xhits().setT(t);
	xhits().setPulse_height(q);
	xhits().setPlane(1);
	xhits().setStrip((int)(10.*x));

	hddm_s::TrdHitList yhits = iter->addTrdHits();
	yhits().setT(t);
	yhits().setPulse_height(q);
	yhits().setPlane(2);
	yhits().setStrip((int)(10.*y));
      }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteTrdTruthHits();
  }
}
