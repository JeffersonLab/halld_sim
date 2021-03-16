#include "CGEMSmearer.h"

//-----------
// cgem_config_t  (constructor)
//-----------
cgem_config_t::cgem_config_t(JEventLoop *loop) 
{
  m_CGEM_WORK_FUNCTION = 20; // 20 pair / eV [FIXME]
  m_CGEM_FANO_FACTOR = 0.1; // [FIXME]
  m_CGEM_ENERGY_THRES = 0.0; // keV 
  m_CGEM_SPATIAL_RESOLUTION = 100; // um
  m_CGEM_TIMING_RESOLUTION = 10; // ns
}

	
//-----------
// SmearEvent
//-----------
void CGEMSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::CgemLayerList layers = record->getCgemLayers();
   hddm_s::CgemLayerList::iterator iter;
   for (iter = layers.begin(); iter != layers.end(); ++iter) {
     
     // If the element already contains a cdcStrawHit list then delete it.
      hddm_s::CgemHitList hits = iter->getCgemHits();
      if (hits.size() > 0) {
         static bool warned = false;
         iter->deleteCgemHits();
         if (!warned) {
	   warned = true;
	   cerr << endl;
	   cerr << "WARNING: CGEM hits already exist in input file! Overwriting!"
		<< endl << endl;
         }
      }
      
     //iter->deleteCGEMHits();
     hddm_s::CgemTruthHitList thits = iter->getCgemTruthHits();
     hddm_s::CgemTruthHitList::iterator titer;
     //layers(0).setLayer(siter->second->layer_);
     for (titer = thits.begin(); titer != thits.end(); ++titer) {
       //int layer = iter->getLayer();
       double E = titer->getDE(); 
       double t = titer->getT(); 
       double x = titer->getX(); 
       double y = titer->getY(); 
       double z = titer->getZ(); 
       //cout <<"Before smearing E " << E << " t " << t << " x " << x << " y " << y << " z " << z << endl; 
       if(config->SMEAR_HITS) {
	 double meanEl = E * 1e9 / cgem_config->m_CGEM_WORK_FUNCTION; // GeV to eV
	 double sigma = sqrt(cgem_config->m_CGEM_FANO_FACTOR * meanEl);
	 int NbEle = (int) (meanEl + gDRandom.SampleGaussian(sigma));
	 //cout << "meanEl " << meanEl << " sigma "  << sigma << " NbEl " << NbEle << endl;
	 double Esmear = (((double) NbEle) * cgem_config->m_CGEM_WORK_FUNCTION) * 1e-6;
	 double tsmear = t + gDRandom.SampleGaussian(cgem_config->m_CGEM_TIMING_RESOLUTION);
	 double xsmear = x + gDRandom.SampleGaussian(cgem_config->m_CGEM_SPATIAL_RESOLUTION * 1e-4);
	 double ysmear = y + gDRandom.SampleGaussian(cgem_config->m_CGEM_SPATIAL_RESOLUTION * 1e-4);
	 double zsmear = z + gDRandom.SampleGaussian(cgem_config->m_CGEM_SPATIAL_RESOLUTION * 1e-4);
	 // Apply a single block threshold. 
	 // Scale threshold by gains
	 if (Esmear >= cgem_config->m_CGEM_ENERGY_THRES){
	   hddm_s::CgemHitList hits = iter->addCgemHits();
	   //hits().setLayer(layer);
	   hits().setDE(Esmear);
	   hits().setT(tsmear);
	   hits().setX(xsmear);
	   hits().setY(ysmear);
	   hits().setZ(zsmear);
	   //cout <<"After smearing E " << Esmear << " t " << tsmear << " x " << xsmear << " y " << ysmear << " z " << zsmear << endl; 
	 }
       }
     }
     
     if (config->DROP_TRUTH_HITS)
       iter->deleteCgemTruthHits();
   }
}

