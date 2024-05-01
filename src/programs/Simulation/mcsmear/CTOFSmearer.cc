#include "CTOFSmearer.h"

//-----------
// ctof_config_t  (constructor)
//-----------
ctof_config_t::ctof_config_t(const std::shared_ptr<const JEvent>& event) 
{
  // default values
  TSIGMA = 100.*k_psec;
  PHOTONS_PERMEV = 400.;
  BAR_THRESHOLD    = 0.0;
  BAR_LENGTH=120.0;
  ATTENUATION_LENGTH=150.;
}


//-----------
// SmearEvent
//-----------
void CTOFSmearer::SmearEvent(hddm_s::HDDM *record)
{
  hddm_s::CtofCounterList ctofs = record->getCtofCounters();
  hddm_s::CtofCounterList::iterator iter;
  for (iter = ctofs.begin(); iter != ctofs.end(); ++iter) {
    // take care of hits
    iter->deleteCtofHits();
    hddm_s::CtofTruthHitList thits = iter->getCtofTruthHits();
    hddm_s::CtofTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      double t = titer->getT();
      double NewE = titer->getDE();
      if(config->SMEAR_HITS) {
	// Smear the time
	t = titer->getT() + gDRandom.SampleGaussian(ctof_config->TSIGMA);
	// Smear the energy
	double npe = titer->getDE() * 1000. * ctof_config->PHOTONS_PERMEV;
	npe += gDRandom.SampleGaussian(sqrt(npe));
	NewE = npe/ctof_config->PHOTONS_PERMEV/1000.;
      }
      // Apply an average attenuation correction to set the energy scale
      NewE *= exp(ctof_config->BAR_LENGTH / 2 / ctof_config->ATTENUATION_LENGTH);
      if (NewE > ctof_config->BAR_THRESHOLD) {
	hddm_s::CtofHitList hits = iter->addCtofHits();
	hits().setEnd(titer->getEnd());
	hits().setT(t);
	hits().setDE(NewE);
      }
    }
    
    if (config->DROP_TRUTH_HITS) {
      iter->deleteCtofTruthHits();
    }
  }
  if (config->DROP_TRUTH_HITS) {
    hddm_s::CppTOFList ctofs = record->getCppTOFs();
    if (ctofs.size() > 0)
      ctofs().deleteCtofTruthPoints();
  }

}
