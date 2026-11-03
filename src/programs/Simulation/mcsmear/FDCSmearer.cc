#include "FDCSmearer.h"

#include <FDC/DFDCGeometry.h>
#include <JANA/JException.h>
#include "DANA/DEvent.h"


#include <sstream>

bool fdc_config_t::FDC_EFFVSDOCA(1);
double fdc_config_t::FDC_CATHODE_SIGMA_SET(-1.);
double fdc_config_t::FDC_DRIFT_SMEAR_NS(0.);
double fdc_config_t::FDC_EFFVSDOCA_PAR_SET[4] = {0.999, 3.75e-4, 0.506, 3.75e-2};

//-----------
// fdc_config_t  (constructor)
//-----------
fdc_config_t::fdc_config_t(const std::shared_ptr<const JEvent>& event) 
{
	// default values
	FDC_TDRIFT_SIGMA      = 0.0;
 	FDC_CATHODE_SIGMA     = 0.0;
 	FDC_PED_NOISE         = 0.0;
 	FDC_THRESHOLD_FACTOR  = 0.0;
 	FDC_TIME_WINDOW       = 0.0;
 	FDC_THRESH_KEV        = 0.0;

	// load data from CCDB
	cout << "Get FDC/fdc_parms parameters from CCDB..." << endl;
    map<string, double> fdcparms;
     if(DEvent::GetCalib(event, "FDC/fdc_parms", fdcparms)) {
     	jerr << "Problem loading FDC/fdc_parms from CCDB!" << endl;
     } else {
       	FDC_TDRIFT_SIGMA      = fdcparms["FDC_TDRIFT_SIGMA"];
       	FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];
       	FDC_THRESHOLD_FACTOR  = fdcparms["FDC_THRESHOLD_FACTOR"];
     	//FDC_PED_NOISE         = fdcparms["FDC_PED_NOISE"];  // ???
       	FDC_TIME_WINDOW       = fdcparms["FDC_TIME_WINDOW"];
       	//FDC_HIT_DROP_FRACTION = fdcparms["FDC_HIT_DROP_FRACTION"];   // ???
       	FDC_THRESH_KEV 		  = fdcparms["FDC_THRESH_KEV"]; 
	}

	// scale for the strip peak amplitude written in fdcDigihit (0 = do not write fdcDigihit)
	FDC_ADC_ASCALE = 0.;
	map<string, double> digi_scales;
	if(DEvent::GetCalib(event, "FDC/digi_scales", digi_scales) || digi_scales.find("FDC_ADC_ASCALE") == digi_scales.end())
		jerr << "Problem loading FDC/digi_scales FDC_ADC_ASCALE from CCDB! fdcDigihit peakAmp will not be written." << endl;
	else
		FDC_ADC_ASCALE = digi_scales["FDC_ADC_ASCALE"];
	jout << "FDC strip peakAmp = q / FDC_ADC_ASCALE, FDC_ADC_ASCALE = " << FDC_ADC_ASCALE << endl;

   	// Calculate ped noise level based on position resolution
   	//   FDC_PED_NOISE = -0.004594 + 0.008711*FDC_CATHODE_SIGMA +
   	//                    0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
	if (FDC_CATHODE_SIGMA_SET >= 0) FDC_CATHODE_SIGMA = FDC_CATHODE_SIGMA_SET;
   	FDC_PED_NOISE = -0.0938 + 0.0485*FDC_CATHODE_SIGMA;
	jout << "FDC_DRIFT_SMEAR_NS = " << FDC_DRIFT_SMEAR_NS << endl;
	jout << "FDC_THRESH_KEV = " << FDC_THRESH_KEV << " (anode wire threshold)" << endl;
	jout << "FDC_CATHODE_SIGMA = " << FDC_CATHODE_SIGMA << " um, FDC_PED_NOISE = " << FDC_PED_NOISE << endl;


	// load efficiency correction factors
	for(int package=1; package<=4; package++) {
		vector< vector<double> > new_strip_efficiencies;
		vector< vector<double> > new_wire_efficiencies;

		char ccdb_str[100];
		sprintf(ccdb_str, "/FDC/package%d/strip_mc_efficiency", package);
    	if(DEvent::GetCalib(event, ccdb_str, new_strip_efficiencies)) {
        	stringstream err_ss;
        	err_ss << "Error loading " << ccdb_str << " !";
        	throw JException(err_ss.str());
        }
		sprintf(ccdb_str,"/FDC/package%d/wire_mc_efficiency", package);
    	if(DEvent::GetCalib(event, ccdb_str, new_wire_efficiencies)) {
        	stringstream err_ss;
        	err_ss << "Error loading " << ccdb_str << " !";
        	throw JException(err_ss.str());
        }

		for(int chamber=0; chamber<6; chamber++) {
        	channel_efficiencies.push_back( new_strip_efficiencies[2*chamber+1] );
        	channel_efficiencies.push_back( new_wire_efficiencies[chamber] );
        	channel_efficiencies.push_back( new_strip_efficiencies[2*chamber] );
		}
	}

    for (int i = 0; i < 4; i++) FDC_EFFVSDOCA_PAR[i] = FDC_EFFVSDOCA_PAR_SET[i];
    jout << "FDC eff vs doca: use " << FDC_EFFVSDOCA << " par " << FDC_EFFVSDOCA_PAR[0] << " " << FDC_EFFVSDOCA_PAR[1]
         << " " << FDC_EFFVSDOCA_PAR[2] << " " << FDC_EFFVSDOCA_PAR[3] << jendl;
}


//-----------
// SmearEvent
//-----------
void FDCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   double t_max = config->TRIGGER_LOOKBACK_TIME + fdc_config->FDC_TIME_WINDOW;
   double threshold = fdc_config->FDC_THRESHOLD_FACTOR * fdc_config->FDC_PED_NOISE; // for sparsification

   hddm_s::FdcChamberList chambers = record->getFdcChambers();
   hddm_s::FdcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {

      // drift-time fluctuation common to the wire hit and the strip hits induced by the same track in this chamber
      // (the strips see the charge induced by the wire avalanche, so the drift delay is the same for both); one shift per itrack
      std::map<int, double> drift_shift_ns;
      auto drift_shift = [&](int itrack) {
         if (!config->SMEAR_HITS || fdc_config->FDC_DRIFT_SMEAR_NS <= 0) return 0.;
         auto it = drift_shift_ns.find(itrack);
         if (it == drift_shift_ns.end())
            it = drift_shift_ns.emplace(itrack, gDRandom.SampleGaussian(fdc_config->FDC_DRIFT_SMEAR_NS)).first;
         return it->second;
      };

      // Add pedestal noise to strip charge data
      hddm_s::FdcCathodeStripList strips = iter->getFdcCathodeStrips();
      hddm_s::FdcCathodeStripList::iterator siter;
      for (siter = strips.begin(); siter != strips.end(); ++siter) {
          // If a fdcCathodeHit already exists delete it
          siter->deleteFdcCathodeHits();
          hddm_s::FdcCathodeTruthHitList thits = 
                                         siter->getFdcCathodeTruthHits();
          hddm_s::FdcCathodeTruthHitList::iterator titer;
          for (titer = thits.begin(); titer != thits.end(); ++titer) {
            // correct simulation efficiencies 
            if (config->APPLY_EFFICIENCY_CORRECTIONS
                  && !gDRandom.DecideToAcceptHit(fdc_config->GetEfficiencyCorrectionFactor(siter)))
              	continue;
          
            double q = titer->getQ();
            double t = titer->getT() + drift_shift(titer->getItrack());
          	if(config->SMEAR_HITS) {
             	q += gDRandom.SampleGaussian(fdc_config->FDC_PED_NOISE);
             	t += gDRandom.SampleGaussian(fdc_config->FDC_TDRIFT_SIGMA)*1.0e9;
			}
            if (q > threshold && t > config->TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcCathodeHitList hits = siter->addFdcCathodeHits();
               hits().setQ(q);
               hits().setT(t);
               // strip peak amplitude on the data scale, read by DEventSourceHDDM as DFDCHit::pulse_height
               if (fdc_config->FDC_ADC_ASCALE > 0) {
                  hddm_s::FdcDigihitList digihit = hits().addFdcDigihits();
                  digihit().setPeakAmp(q / fdc_config->FDC_ADC_ASCALE);
               }
            }
         }

         if (config->DROP_TRUTH_HITS)
            siter->deleteFdcCathodeTruthHits();
      }

      // Add drift time varation to the anode data 
      hddm_s::FdcAnodeWireList wires = iter->getFdcAnodeWires();
      hddm_s::FdcAnodeWireList::iterator witer;
      for (witer = wires.begin(); witer != wires.end(); ++witer) {
         // If a fdcAnodeHit exists already delete it
         witer->deleteFdcAnodeHits();
         hddm_s::FdcAnodeTruthHitList thits = witer->getFdcAnodeTruthHits();
         hddm_s::FdcAnodeTruthHitList::iterator titer;
         for (titer = thits.begin(); titer != thits.end(); ++titer) {
            // anode discriminator threshold, applied to the energy deposited on the wire
            // (dE in GeV, summed by hdgeant4 over the track segments it combined into this hit)
            if (titer->getDE() * 1e6 < fdc_config->FDC_THRESH_KEV)
               continue;
             // correct simulation efficiencies 
		     if (config->APPLY_EFFICIENCY_CORRECTIONS
             		&& !gDRandom.DecideToAcceptHit(fdc_config->GetEfficiencyCorrectionFactor(witer)))
             		continue;
            double doca = titer->getD();
            if (config->APPLY_EFFICIENCY_CORRECTIONS
                    && !gDRandom.DecideToAcceptHit(fdc_config->GetEfficiencyVsDOCA(doca)))
               continue;

            double t = titer->getT() + drift_shift(titer->getItrack());
          	if(config->SMEAR_HITS) {
               t += gDRandom.SampleGaussian(fdc_config->FDC_TDRIFT_SIGMA)*1.0e9;
            }
            if (t > config->TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcAnodeHitList hits = witer->addFdcAnodeHits();
               hits().setT(t);
               hits().setDE(titer->getDE());
            }
         }

         if (config->DROP_TRUTH_HITS)
            witer->deleteFdcAnodeTruthHits();
      }
      if (config->DROP_TRUTH_HITS)
         iter->deleteFdcTruthPoints();
   }
}

