#include "ITOFSmearer.h"

//-----------
// itof_config_t  (constructor)
//-----------
itof_config_t::itof_config_t(JEventLoop *loop) 
{
  // default values
  TSIGMA = 250.*k_psec;
  PHOTONS_PERMEV = 400.;
}


//-----------
// SmearEvent
//-----------
void ITOFSmearer::SmearEvent(hddm_s::HDDM *record)
{
  hddm_s::ItofCounterList itofs = record->getItofCounters();
  hddm_s::ItofCounterList::iterator iter;
  for (iter = itofs.begin(); iter != itofs.end(); ++iter) {
    // take care of hits
    iter->deleteItofHits();
    hddm_s::ItofTruthHitList thits = iter->getItofTruthHits();
    hddm_s::ItofTruthHitList::iterator titer;
    for (titer = thits.begin(); titer != thits.end(); ++titer) {
      double t = titer->getT();
      double NewE = titer->getDE();
      if(config->SMEAR_HITS) {
	// Smear the time
	t = titer->getT() + gDRandom.SampleGaussian(itof_config->TSIGMA);
	// Smear the energy
	double npe = titer->getDE() * 1000. * itof_config->PHOTONS_PERMEV;
	npe += gDRandom.SampleGaussian(sqrt(npe));
	NewE = npe/itof_config->PHOTONS_PERMEV/1000.;
      }
      {
      hddm_s::ItofHitList hits = iter->addItofHits();
	hits().setT(t);
	hits().setDE(NewE);
	hits().setX(titer->getX());
	hits().setY(titer->getY());
      }
    }
    
    if (config->DROP_TRUTH_HITS) {
      iter->deleteItofTruthHits();
    }
  }
  if (config->DROP_TRUTH_HITS) {
    hddm_s::InnerTOFList itofs = record->getInnerTOFs();
    if (itofs.size() > 0)
      itofs().deleteItofTruthPoints();
  }

}
