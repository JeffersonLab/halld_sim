#include "GEMTRDSmearer.h"

//-----------
// gemtrd_config_t  (constructor)
//-----------
gemtrd_config_t::gemtrd_config_t(const std::shared_ptr<const JEvent>& event) 
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

    // Keep lists of strips that have hits
    vector<strip_hit_t>x_strip_hits;
    vector<strip_hit_t>y_strip_hits;
    
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
      }
      
      // Distribute charge over strips
      int strip_y0=528/2-int(10.*(y+59.031));
      int strip_x0=720/2-int(10*(x+47.4695));
      for (int strip=-10;strip<=10;strip++){
	// Assume a symmetric charge distribution in x and y
	double q_strip=GetStripCharge(q,strip);
	
	int mystrip=strip+strip_x0;
	if (mystrip>0&&mystrip<721){
	  x_strip_hits.push_back(strip_hit_t(1,mystrip,q_strip,t));
	}
	mystrip=strip+strip_y0;
	if (mystrip>0&&mystrip<433){
	  y_strip_hits.push_back(strip_hit_t(2,mystrip,q_strip,t));
	}
      }
    }
    vector<bool>x_used(x_strip_hits.size());
    for (unsigned int i=0;i<x_strip_hits.size();i++){
      if (!x_used[i]){
	strip_hit_t myhit=x_strip_hits[i];
	for (unsigned int j=i+1;j<x_strip_hits.size();j++){
	  strip_hit_t myhit2=x_strip_hits[j];
	  if (myhit.strip==myhit2.strip && fabs(myhit.t-myhit2.t)<8.){
	    myhit.q+=myhit2.q;
	    if (myhit2.t<myhit.t) myhit.t=myhit2.t;
	    x_used[j]=true;
	  }
	}
      }      
    }
    for (unsigned int i=0;i<x_strip_hits.size();i++){
      strip_hit_t myhit=x_strip_hits[i];
      if (x_used[i]==false&&myhit.q>0.01){
	hddm_s::GemtrdHitList hits = iter->addGemtrdHits();
	double t=8.0*floor(myhit.t/8.0);
	hits().setT(t);
	hits().setQ(myhit.q);
	hits().setPlane(1);
	hits().setStrip(myhit.strip);
      }
    }
    vector<bool>y_used(y_strip_hits.size());
    for (unsigned int i=0;i<y_strip_hits.size();i++){
      if (!y_used[i]){
	strip_hit_t myhit=y_strip_hits[i];
	for (unsigned int j=i+1;j<y_strip_hits.size();j++){
	  strip_hit_t myhit2=y_strip_hits[j];
	  if (myhit.strip==myhit2.strip && fabs(myhit.t-myhit2.t)<8.){
	    myhit.q+=myhit2.q;
	    if (myhit2.t<myhit.t) myhit.t=myhit2.t;
	    y_used[j]=true;
	  }
	}
      }      
    }
    for (unsigned int i=0;i<y_strip_hits.size();i++){
      strip_hit_t myhit=y_strip_hits[i];
      if (y_used[i]==false && myhit.q>0.01){
	hddm_s::GemtrdHitList hits = iter->addGemtrdHits();
	double t=8.0*floor(myhit.t/8.0);
	hits().setT(t);
	hits().setQ(myhit.q);
	hits().setPlane(2);
	hits().setStrip(myhit.strip);
      }
    }
    
    if (config->DROP_TRUTH_HITS)
      iter->deleteGemtrdTruthHits();
  }
}

// Induce charge on the strips according to the semi-empirical Mathieson
// function
double GEMTRDSmearer::GetStripCharge(double q,int strip) const {
  const double K2=5.0; // induced charge distribution shape parameter
  double factor= 0.25 * M_PI * K2;
  double strip_pitch=0.1; // cm
  double strip_gap=0.01; //cm, place holder
  double anode_cathode_separation=0.5; // cm
  double lambda1 = ((strip - 0.5) * strip_pitch +
		    0.5*strip_gap) / anode_cathode_separation;
  double lambda2 = ((strip + 0.5) * strip_pitch - 
		    0.5*strip_gap) / anode_cathode_separation;
  return 0.25 * q * (tanh(factor * lambda2) - tanh(factor * lambda1));
}
