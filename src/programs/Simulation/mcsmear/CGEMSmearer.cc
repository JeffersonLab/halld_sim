#include "CGEMSmearer.h"

//-----------
// cgem_config_t  (constructor)
//-----------
cgem_config_t::cgem_config_t(JEventLoop *loop) 
{
                  
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
       
       //if(config->SMEAR_HITS) {
       //	}
       
       // Apply a single block threshold. 
       // Scale threshold by gains
       //if (E >= Ethreshold){
       hddm_s::CgemHitList hits = iter->addCgemHits();
       //hits().setLayer(layer);
       hits().setDE(E);
       hits().setT(t);
       hits().setX(x);
       hits().setY(y);
       hits().setZ(z);
     }
      
     if (config->DROP_TRUTH_HITS)
       iter->deleteCgemTruthHits();
   }
}
