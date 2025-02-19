#include "TOFSmearer.h"

#include <TOF/DTOFGeometry.h>
#include "DANA/DEvent.h"



//-----------
// tof_config_t  (constructor)
//-----------
tof_config_t::tof_config_t(const std::shared_ptr<const JEvent>& event) 
{
	// default values
 	TOF_SIGMA = 100.*k_psec;
 	TOF_PHOTONS_PERMEV = 400.;
 	TOF_BAR_THRESHOLD    = 0.0;
	ATTENUATION_LENGTH = 150.;
	FULL_BAR_LENGTH = 252.;

	// load values from geometry
	vector <const DTOFGeometry*> TOFGeom;
	event->Get(TOFGeom);
	TOF_NUM_PLANES = TOFGeom[0]->Get_NPlanes();
	TOF_NUM_BARS = TOFGeom[0]->Get_NBars();

	// Load data from CCDB
    string locTOFParmsTable = TOFGeom[0]->Get_CCDB_DirectoryName() + "/tof_parms";
    cout<<"Get "<<locTOFParmsTable<<" parameters from CCDB..."<<endl;
    map<string, double> tofparms;
    if(DEvent::GetCalib(event, locTOFParmsTable.c_str(), tofparms)) {
     	jerr << "Problem loading "<<locTOFParmsTable<<" from CCDB!" << endl;
     	return;
    }
     	
    TOF_SIGMA =  tofparms.at("TOF_SIGMA");
    TOF_PHOTONS_PERMEV =  tofparms.at("TOF_PHOTONS_PERMEV");
    ATTENUATION_LENGTH = tofparms.at("TOF_ATTEN_LENGTH");
    FULL_BAR_LENGTH = tofparms.at("TOF_PADDLE_LENGTH");
	
    string locTOFPaddleResolTable = TOFGeom[0]->Get_CCDB_DirectoryName() + "/paddle_resolutions";
	cout<<"get "<<locTOFPaddleResolTable<<" from calibDB"<<endl;
    vector <double> TOF_PADDLE_TIME_RESOLUTIONS_TEMP;
    if(DEvent::GetCalib(event, locTOFPaddleResolTable.c_str(), TOF_PADDLE_TIME_RESOLUTIONS_TEMP)) {
    	jerr << "Problem loading "<<locTOFPaddleResolTable<<" from CCDB!" << endl;
    } else {
    	for (unsigned int i = 0; i < TOF_PADDLE_TIME_RESOLUTIONS_TEMP.size(); i++) {
       		TOF_PADDLE_TIME_RESOLUTIONS.push_back(TOF_PADDLE_TIME_RESOLUTIONS_TEMP.at(i));
    	}
    }

	// load per-channel efficiencies
    string locTOFChannelEffTable = TOFGeom[0]->Get_CCDB_DirectoryName() + "/channel_mc_efficiency";
	vector<double> raw_table;
	if(DEvent::GetCalib(event, locTOFChannelEffTable.c_str(), raw_table)) {
    	jerr << "Problem loading "<<locTOFChannelEffTable<<" from CCDB!" << endl;
    } else {
    	int channel = 0;

    	for(int plane=0; plane<TOF_NUM_PLANES; plane++) {
        	int plane_index=2*TOF_NUM_BARS*plane;
        	channel_efficiencies.push_back( vector< pair<double,double> >(TOF_NUM_BARS) );
        	for(int bar=0; bar<TOF_NUM_BARS; bar++) {
            	channel_efficiencies[plane][bar] 
                	= pair<double,double>(raw_table[plane_index+bar],
                    	    raw_table[plane_index+TOF_NUM_BARS+bar]);
            	channel+=2;
            }	      
        }
    }

	cout << "Number of TOF bars per plane = " << TOF_NUM_BARS << endl;
}


//-----------
// SmearEvent
//-----------
void TOFSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::FtofCounterList tofs = record->getFtofCounters();
   hddm_s::FtofCounterList::iterator iter;
   for (iter = tofs.begin(); iter != tofs.end(); ++iter) {
      // take care of hits
      iter->deleteFtofHits();
      hddm_s::FtofTruthHitList thits = iter->getFtofTruthHits();
      hddm_s::FtofTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
         // correct simulation efficiencies 
		 if (config->APPLY_EFFICIENCY_CORRECTIONS
		 		&& !gDRandom.DecideToAcceptHit(tof_config->GetEfficiencyCorrectionFactor(titer)))
		 			continue;
		 			
         // Smear the time
         //double t = titer->getT() + gDRandom.SampleGaussian(tof_config->TOF_SIGMA);
         double t = titer->getT();
         // Smear the energy
         float NewE = titer->getDE();
         if(config->SMEAR_HITS) {
			 t += gDRandom.SampleGaussian(tof_config->GetHitTimeResolution(iter->getPlane(),iter->getBar()));
         	 double npe = titer->getDE() * 1000. * tof_config->TOF_PHOTONS_PERMEV;
         	 npe += gDRandom.SampleGaussian(sqrt(npe));
          	 NewE = npe/tof_config->TOF_PHOTONS_PERMEV/1000.;
		 }
         // Apply an average attenuation correction to set the energy scale
         NewE *= exp(tof_config->FULL_BAR_LENGTH / 2 / tof_config->ATTENUATION_LENGTH);
         if (NewE > tof_config->TOF_BAR_THRESHOLD) {
            hddm_s::FtofHitList hits = iter->addFtofHits();
            hits().setEnd(titer->getEnd());
            hits().setT(t);
            hits().setDE(NewE);
         }
      }
    
      if (config->DROP_TRUTH_HITS) {
         iter->deleteFtofTruthHits();
      }
   }
   if (config->DROP_TRUTH_HITS) {
      hddm_s::ForwardTOFList tofs = record->getForwardTOFs();
      if (tofs.size() > 0)
         tofs().deleteFtofTruthPoints();
   }
}
