#include "GEMTRDSmearer.h"

//-----------
// gemtrd_config_t  (constructor)
//-----------
gemtrd_config_t::gemtrd_config_t(const std::shared_ptr<const JEvent>& event) 
{
  // default values
  GEMTRD_TSIGMA = 1.0;  // ns
  GEMTRD_ASIGMA = 1.0; // ???
  GEMTRD_XYSIGMA = 0.01; // cm
  GEMTRD_THRESHOLD = 10.;
  GEMTRD_INTEGRAL_TO_AMPLITUDE=1./28.8; // copied from CDC
  GEMTRD_INSTALLED=true;
  
  map<string,string> installed;
  DEvent::GetCalib(event, "/TRD/install_status", installed);
  if(atoi(installed["status"].data()) == 0){
    GEMTRD_INSTALLED=false;
  }
  else {
    // Get the geometry
    DGeometry* locGeometry = DEvent::GetDGeometry(event);
    locGeometry->GetGEMTRDxy_vec(GEMTRDx,GEMTRDy);
  }
}

//-----------
// SmearEvent
//-----------
void GEMTRDSmearer::SmearEvent(hddm_s::HDDM *record)
{
  if (gemtrd_config->GEMTRD_INSTALLED==false) return;
      
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
    int chamber=iter->getLayer()-1;
    
    // Keep lists of strips that have hits
    map<int,vector<pair<double,double>>>x_strip_hits;
    map<int,vector<pair<double,double>>>y_strip_hits;
    
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
      int strip_y0=528/2-int(10.*(y-gemtrd_config->GEMTRDy[chamber]));
      if (chamber==1){
	strip_y0=528/2+int(10.*(y-gemtrd_config->GEMTRDy[chamber]));
      }
      int strip_x0=720/2-int(10*(x-gemtrd_config->GEMTRDx[chamber]));
      for (int strip=-10;strip<=10;strip++){
	// Assume a symmetric charge distribution in x and y
	double q_strip=GetStripCharge(q,strip);
	
	int mystrip=strip+strip_x0;
	if (mystrip>0&&mystrip<721){
	  x_strip_hits[mystrip].push_back(make_pair(q_strip,t));
	}
	mystrip=strip+strip_y0;
	if (mystrip>0&&mystrip<529){
	  y_strip_hits[mystrip].push_back(make_pair(q_strip,t));
	}
      }
    }
    
    for (auto siter=x_strip_hits.begin(); siter!=x_strip_hits.end(); ++siter){
      MakeHits(1,siter->first,siter->second,iter);
    }
    for (auto siter=y_strip_hits.begin(); siter!=y_strip_hits.end(); ++siter){
      MakeHits(2,siter->first,siter->second,iter);
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

double GEMTRDSmearer::EmulateSignal(double t,vector<pair<double,double>>&hits) const{
   double asic_gain = 2.3; // mv/fC
   double signal_mV = 0;
   for (auto iter = hits.begin(); iter != hits.end(); ++iter){
     
     if (t>iter->second){
       double my_dt=t-iter->second;
       signal_mV+=asic_gain * iter->first * AsicResponse(my_dt);
     }
   }
   return signal_mV;
}
double GEMTRDSmearer::AsicResponse(double t_ns) const {
   // Simulation of the ASIC response to a pulse due to a cluster
  double par[11] = {-0.01986, 0.01802, -0.001097, 10.3, 11.72,
                     -0.03701, 35.84, 15.93, 0.006141, 80.95, 24.77};
  if (t_ns < par[3])
    return par[0] * t_ns + par[1] * t_ns*t_ns + par[2] * t_ns*t_ns*t_ns;
  else
    return (par[0] * par[3] +
	    par[1] * par[3]*par[3] +
	    par[2] * par[3]*par[3]*par[3]) *
      exp(-pow((t_ns - par[3])/par[4], 2)) +
      par[5] * exp(-pow((t_ns - par[6])/par[7], 2)) +
      par[8] * exp(-pow((t_ns - par[9])/par[10], 2));
}

void GEMTRDSmearer::MakeHits(int plane,int strip,
			     vector<pair<double,double>>&data,
			     hddm_s::GemtrdChamberList::iterator iter
			     ) const {
  // Emulate pulse signal
  vector<int>nanosamples(1000);
  for (int i=0;i<1000;i++){
    nanosamples[i]=EmulateSignal(double(i),data);
  }
  vector<double>samples(125);
  for (int i=0;i<1000;i+=8){
    double sum=0;
    for (int j=i;j<i+8;j++){
      sum+=nanosamples[j];
    }
    sum+=gDRandom.SampleGaussian(gemtrd_config->GEMTRD_ASIGMA);
    samples[i]=sum;
  }
  // Emulate FADC sampling
  int num_peaks=0;
  for (int i=1;i<124;i++){
    if (samples[i]>gemtrd_config->GEMTRD_THRESHOLD){
      if (samples[i-1]<samples[i] && samples[i+1]<samples[i]){
	// found peak
	hddm_s::GemtrdHitList hits = iter->addGemtrdHits();
	double t=8.*double(i);
	hits().setChamber(iter->getLayer());
	hits().setT(t);
	hits().setQ(samples[i]);
	hits().setPlane(plane);
	hits().setStrip(strip);

	num_peaks++;
      }
    }
    if (num_peaks==15) break;
  }
}
