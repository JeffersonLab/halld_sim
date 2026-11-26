#include "FDCSmearer.h"

#include <FDC/DFDCGeometry.h>
#include <JANA/JException.h>
#include "DANA/DEvent.h"


#include <sstream>
#include <fstream>

bool fdc_config_t::FDC_EFFVSDOCA(1);
double fdc_config_t::FDC_CATHODE_SIGMA_SET(-1.);
double fdc_config_t::FDC_DRIFT_SMEAR_NS(0.);
double fdc_config_t::FDC_GAIN_SCALE(1.);
double fdc_config_t::FDC_GAIN_PKG[4] = {1., 1., 1., 1.};
double fdc_config_t::FDC_DOCA_ATTEN(0.);
double fdc_config_t::FDC_GAIN_SIGMA(0.);
double fdc_config_t::FDC_STRIP_THR_SET(-1.);
double fdc_config_t::FDC_WIRE_THR_KEV_SET(-1.);
double fdc_config_t::FDC_EFFVSDOCA_PAR_SET[4] = {0.999, 3.75e-4, 0.506, 3.75e-2};
double fdc_config_t::FDC_RAW_AMP_THR(-1.);
const double fdc_config_t::FDC_RAW_AMP_THR_DEFAULT(50.);
double fdc_config_t::FDC_PEAKAMP_ASCALE(-1.);
std::string fdc_config_t::FDC_Q4_TABLE("");
std::string fdc_config_t::FDC_STRIP_TIME_TABLE("");
double fdc_config_t::FDC_STRIP_T0_PKG[4] = {0., 0., 0., 0.};
double fdc_config_t::FDC_REST_ASCALE(-1.);

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
	if (FDC_PEAKAMP_ASCALE > 0)
		jout << "FDC strip peakAmp = q / FDC:PEAKAMP_ASCALE = q / " << FDC_PEAKAMP_ASCALE
		     << " in every run; readout threshold on q / FDC_ADC_ASCALE = q / " << FDC_ADC_ASCALE << endl;

   	// Calculate ped noise level based on position resolution
   	//   FDC_PED_NOISE = -0.004594 + 0.008711*FDC_CATHODE_SIGMA +
   	//                    0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
	if (FDC_CATHODE_SIGMA_SET >= 0) FDC_CATHODE_SIGMA = FDC_CATHODE_SIGMA_SET;
   	FDC_PED_NOISE = -0.0938 + 0.0485*FDC_CATHODE_SIGMA;
	jout << "FDC_DRIFT_SMEAR_NS = " << FDC_DRIFT_SMEAR_NS << endl;
	jout << "FDC_GAIN_SCALE = " << FDC_GAIN_SCALE << " x package " << FDC_GAIN_PKG[0] << "/" << FDC_GAIN_PKG[1] << "/" << FDC_GAIN_PKG[2] << "/" << FDC_GAIN_PKG[3] << ", FDC_DOCA_ATTEN = " << FDC_DOCA_ATTEN << ", FDC_GAIN_SIGMA = " << FDC_GAIN_SIGMA << ", FDC_STRIP_THRESHOLD = " << FDC_STRIP_THR_SET << " pC, FDC_WIRE_THRESHOLD_KEV = " << FDC_WIRE_THR_KEV_SET << " (< 0: CCDB)" << endl;
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

	// readout thresholds of the cathode strips on raw_amp = q / FDC_ADC_ASCALE, in ADC counts, per strip from
	// /FDC/packageN/strip_thresholds (12 rows, as strip_mc_efficiency); FDC_RAW_AMP_THR_DEFAULT where a table is absent
	for(int package=1; package<=4; package++) {
		vector< vector<double> > new_thresholds;
		char ccdb_str[100];
		sprintf(ccdb_str, "/FDC/package%d/strip_thresholds", package);
		bool have_table = false;
		if (FDC_RAW_AMP_THR < 0) {
			have_table = !DEvent::GetCalib(event, ccdb_str, new_thresholds) && new_thresholds.size() == 12;
			if (have_table)
				jout << "FDC strip raw-amplitude thresholds of package " << package << " from " << ccdb_str << endl;
			else
				jout << "FDC strip raw-amplitude thresholds of package " << package << ": no " << ccdb_str
				     << ", default " << FDC_RAW_AMP_THR_DEFAULT << " ADC counts" << endl;
		}
		for(int chamber=0; chamber<6; chamber++) {
			channel_raw_amp_thresholds.push_back( have_table? new_thresholds[2*chamber+1] : vector<double>() );
			channel_raw_amp_thresholds.push_back( vector<double>() );
			channel_raw_amp_thresholds.push_back( have_table? new_thresholds[2*chamber] : vector<double>() );
		}
	}
	if (FDC_RAW_AMP_THR >= 0)
		jout << "FDC strip raw-amplitude threshold FDC:RAW_AMP_THRESHOLD = " << FDC_RAW_AMP_THR << " ADC counts for every strip" << endl;

	// q4 per-channel table (FDC:Q4_TABLE); every strip, wire and wire plane must be given, else the run is refused
	q4 = !FDC_Q4_TABLE.empty();
	if (q4) {
		for (int i = 0; i < 48; i++) for (int j = 0; j < 192; j++) q4_gain[i][j] = q4_thr[i][j] = -1;
		for (int i = 0; i < 24; i++) { q4_q50[i] = q4_width[i] = -1; for (int j = 0; j < 96; j++) q4_eff[i][j] = -1; }
		std::ifstream tf(FDC_Q4_TABLE.c_str());
		if (!tf) throw JException("FDC:Q4_TABLE " + FDC_Q4_TABLE + " cannot be opened");
		std::string line; int ns = 0, nw = 0, nt = 0;
		while (std::getline(tf, line)) {
			std::istringstream ls(line); std::string k; if (!(ls >> k) || k[0] == '#') continue;
			int a, b; double x, y;
			if (k == "S" && (ls >> a >> b >> x >> y) && a >= 1 && a <= 48 && b >= 1 && b <= 192 && x > 0) { q4_gain[a-1][b-1] = x; q4_thr[a-1][b-1] = y; ns++; }
			else if (k == "W" && (ls >> a >> b >> x) && a >= 1 && a <= 24 && b >= 1 && b <= 96 && x >= 0 && x <= 1) { q4_eff[a-1][b-1] = x; nw++; }
			else if (k == "T" && (ls >> a >> x >> y) && a >= 1 && a <= 24 && y > 0) { q4_q50[a-1] = x; q4_width[a-1] = y; nt++; }
			else throw JException("FDC:Q4_TABLE " + FDC_Q4_TABLE + ": bad line: " + line);
		}
		for (int i = 0; i < 48; i++) for (int j = 0; j < 192; j++) if (q4_gain[i][j] <= 0) throw JException("FDC:Q4_TABLE " + FDC_Q4_TABLE + ": strip missing");
		for (int i = 0; i < 24; i++) { if (q4_width[i] <= 0) throw JException("FDC:Q4_TABLE " + FDC_Q4_TABLE + ": wire turn-on missing");
			for (int j = 0; j < 96; j++) if (q4_eff[i][j] < 0) throw JException("FDC:Q4_TABLE " + FDC_Q4_TABLE + ": wire missing"); }
		double gs = 0, ts = 0; for (int i = 0; i < 48; i++) for (int j = 0; j < 192; j++) { gs += q4_gain[i][j]; ts += q4_thr[i][j]; }
		jout << "FDC q4 table " << FDC_Q4_TABLE << ": " << ns << " strips (mean gain " << gs / 9216 << " ADC/pC, mean threshold " << ts / 9216
		     << " ADC), " << nw << " wires, " << nt << " wire turn-on planes" << jendl;
	}
	// q5 strip time table (FDC:STRIP_TIME_TABLE)
	if (!FDC_STRIP_TIME_TABLE.empty()) {
		std::ifstream sf(FDC_STRIP_TIME_TABLE.c_str());
		if (!sf) throw JException("FDC:STRIP_TIME_TABLE " + FDC_STRIP_TIME_TABLE + " cannot be opened");
		std::string line;
		while (std::getline(sf, line)) {
			if (line.empty() || line[0] == '#') continue;
			if (line[0] == 'Q') {   // quantile form: Q <centre> <q_0 .. q_N>
				std::istringstream ls(line.substr(1)); double ce, v; std::vector<double> q;
				if (!(ls >> ce) || ce <= 0 || (!stt_lnA.empty() && log(ce) <= stt_lnA.back()))
					throw JException("FDC:STRIP_TIME_TABLE " + FDC_STRIP_TIME_TABLE + ": bad line: " + line.substr(0, 60));
				while (ls >> v) q.push_back(std::max(v, -100.));   // the lowest quantile is the histogram edge
				if (q.size() < 11 || (!stt_q.empty() && q.size() != stt_q[0].size()))
					throw JException("FDC:STRIP_TIME_TABLE " + FDC_STRIP_TIME_TABLE + ": bad quantile count: " + line.substr(0, 60));
				for (size_t k = 1; k < q.size(); k++) if (q[k] < q[k-1]) throw JException("FDC:STRIP_TIME_TABLE: quantiles not ordered");
				stt_lnA.push_back(log(ce)); stt_q.push_back(q); stt_shift.push_back(q[q.size() / 2]); stt_width.push_back(0.);
				continue;
			}
			if (!stt_q.empty()) throw JException("FDC:STRIP_TIME_TABLE: mixed table forms");
			std::istringstream ls(line); double lo, hi, ce, sh, w;
			if (!(ls >> lo >> hi >> ce >> sh >> w) || ce <= 0 || w < 0 || (!stt_lnA.empty() && log(ce) <= stt_lnA.back()))
				throw JException("FDC:STRIP_TIME_TABLE " + FDC_STRIP_TIME_TABLE + ": bad line: " + line);
			stt_lnA.push_back(log(ce)); stt_shift.push_back(sh); stt_width.push_back(w);
		}
		if (stt_lnA.empty()) throw JException("FDC:STRIP_TIME_TABLE " + FDC_STRIP_TIME_TABLE + ": no rows");
		if (!q4) throw JException("FDC:STRIP_TIME_TABLE needs FDC:Q4_TABLE (the strip ADC)");
		if (!stt_q.empty())
			jout << "FDC strip time quantile table " << FDC_STRIP_TIME_TABLE << ": " << stt_lnA.size() << " rows x " << stt_q[0].size()
			     << " quantiles, median " << stt_shift.front() << " ns at " << exp(stt_lnA.front()) << " ADC to " << stt_shift.back() << " ns at "
			     << exp(stt_lnA.back()) << " ADC" << jendl;
		else
		jout << "FDC strip time table " << FDC_STRIP_TIME_TABLE << ": " << stt_lnA.size() << " rows, shift " << stt_shift.front() << " ns at "
		     << exp(stt_lnA.front()) << " ADC to " << stt_shift.back() << " ns at " << exp(stt_lnA.back()) << " ADC, width " << stt_width.front()
		     << " to " << stt_width.back() << " ns" << jendl;
	}
	if (FDC_STRIP_T0_PKG[0] != 0 || FDC_STRIP_T0_PKG[1] != 0 || FDC_STRIP_T0_PKG[2] != 0 || FDC_STRIP_T0_PKG[3] != 0) {
		if (!q4) throw JException("FDC:STRIP_T0_PKG1..4 need FDC:Q4_TABLE");
		jout << "FDC strip time residue per package: " << FDC_STRIP_T0_PKG[0] << " " << FDC_STRIP_T0_PKG[1] << " " << FDC_STRIP_T0_PKG[2]
		     << " " << FDC_STRIP_T0_PKG[3] << " ns (subtracted from the strip times)" << jendl;
	}
	if (FDC_REST_ASCALE > 0)
		jout << "FDC:REST_ASCALE = " << FDC_REST_ASCALE << ": strip q = this x peakAmp after the random-trigger merge" << jendl;

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
   if (fdc_config->FDC_STRIP_THR_SET >= 0) threshold = fdc_config->FDC_STRIP_THR_SET;
   double wire_thr_kev = (fdc_config->FDC_WIRE_THR_KEV_SET >= 0) ? fdc_config->FDC_WIRE_THR_KEV_SET : fdc_config->FDC_THRESH_KEV;
   double gain_global = fdc_config->FDC_GAIN_SCALE;   // gas gain relative to the simulation: scales the avalanche charge seen on wires and strips

   hddm_s::FdcChamberList chambers = record->getFdcChambers();
   hddm_s::FdcChamberList::iterator iter;
   for (iter = chambers.begin(); iter != chambers.end(); ++iter) {
      int pkg = (iter->getModule() - 1) / 2;   // FDC package 0-3 (3 planes per module, 2 modules per package)
      double gain = gain_global * ((pkg >= 0 && pkg < 4) ? fdc_config->FDC_GAIN_PKG[pkg] : 1.);
      const bool q4 = fdc_config->q4;
      if (q4) gain = 1.;   // q4: the gain is per strip, in the table (ADC per pC)
      // q4: strip ADC summed per track on each cathode of this chamber ([0] plane 1, [1] plane 3), for the wire turn-on
      std::map<int, std::pair<double, double> > q4_qsum;

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

      // signal attenuation vs the track's doca to its wire in this chamber (charge lost at the cell edge), per itrack, from the anode truth hits
      std::map<int, double> atten;
      auto atten_of_doca = [&](double d) { return std::max(0.05, 1. - fdc_config->FDC_DOCA_ATTEN * d * d); };
      if (fdc_config->FDC_DOCA_ATTEN > 0) {
         hddm_s::FdcAnodeWireList awires = iter->getFdcAnodeWires();
         for (auto aw = awires.begin(); aw != awires.end(); ++aw) {
            hddm_s::FdcAnodeTruthHitList ath = aw->getFdcAnodeTruthHits();
            for (auto at = ath.begin(); at != ath.end(); ++at) {
               double a = atten_of_doca(at->getD());
               auto it = atten.find(at->getItrack());
               if (it == atten.end() || a > it->second) atten[at->getItrack()] = a;   // the wire with the closest approach
            }
         }
      }
      auto atten_of_track = [&](int itrack) { auto it = atten.find(itrack); return (it == atten.end()) ? 1. : it->second; };
      // avalanche gain fluctuation, one factor per track in this chamber (the strips see the charge induced by the same wire avalanche)
      std::map<int, double> fluct;
      auto fluct_of_track = [&](int itrack) {
         if (fdc_config->FDC_GAIN_SIGMA <= 0) return 1.;
         auto it = fluct.find(itrack);
         if (it == fluct.end()) {
            double s = fdc_config->FDC_GAIN_SIGMA;
            it = fluct.emplace(itrack, exp(gDRandom.SampleGaussian(s) - 0.5 * s * s)).first;
         }
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
          
            double q = titer->getQ() * gain * atten_of_track(titer->getItrack()) * fluct_of_track(titer->getItrack());
            double t = titer->getT() + drift_shift(titer->getItrack());
          	if(config->SMEAR_HITS) {
             	q += gDRandom.SampleGaussian(fdc_config->FDC_PED_NOISE);
             	t += gDRandom.SampleGaussian(fdc_config->FDC_TDRIFT_SIGMA)*1.0e9;
			}
            if (q4) {
               // q4: strip ADC = gain x q, readout threshold per strip, both from the table
               int cp = fdc_config_t::q4_cplane(siter->getModule(), siter->getLayer(), siter->getPlane());
               int s = siter->getStrip(); if (s > 1000) s -= 1000;
               if (cp < 1 || cp > 48 || s < 1 || s > 192) continue;
               double adc = fdc_config->q4_gain[cp-1][s-1] * q;
               if (config->SMEAR_HITS && !fdc_config->stt_lnA.empty()) {
                  // q5: strip time walk and spread vs strip amplitude, measured from data (FDC:STRIP_TIME_TABLE)
                  if (!fdc_config->stt_q.empty()) t += fdc_config->StripTimeDraw(adc, gDRandom.Rndm());
                  else { double sh, w; fdc_config->StripTime(adc, sh, w); t += sh + gDRandom.SampleGaussian(w); }
               }
               t -= fdc_config_t::FDC_STRIP_T0_PKG[(cp - 1) / 12];   // q5 item 5: strip time residue of the data REST, per package
               if (adc < fdc_config->q4_thr[cp-1][s-1] || t <= config->TRIGGER_LOOKBACK_TIME || t >= t_max)
                  continue;
               hddm_s::FdcCathodeHitList hits = siter->addFdcCathodeHits();
               hits().setQ(q);
               hits().setT(t);
               hddm_s::FdcDigihitList digihit = hits().addFdcDigihits();
               digihit().setPeakAmp(adc);
               auto &qs = q4_qsum[titer->getItrack()];
               if (siter->getPlane() == 3) qs.second += adc; else qs.first += adc;
               continue;
            }
            // readout threshold on the raw amplitude raw_amp = q / FDC_ADC_ASCALE, in ADC counts
            if (fdc_config->FDC_ADC_ASCALE > 0
                  && q / fdc_config->FDC_ADC_ASCALE < fdc_config->GetRawAmpThreshold(siter))
               continue;
            if (q > threshold && t > config->TRIGGER_LOOKBACK_TIME && t < t_max) {
               hddm_s::FdcCathodeHitList hits = siter->addFdcCathodeHits();
               hits().setQ(q);
               hits().setT(t);
               // strip peak amplitude on the data scale, read by DEventSourceHDDM as DFDCHit::pulse_height
               if (fdc_config->FDC_ADC_ASCALE > 0) {
                  hddm_s::FdcDigihitList digihit = hits().addFdcDigihits();
                  digihit().setPeakAmp(q / ((fdc_config->FDC_PEAKAMP_ASCALE > 0)? fdc_config->FDC_PEAKAMP_ASCALE
                                                                                : fdc_config->FDC_ADC_ASCALE));
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
            double att = (fdc_config->FDC_DOCA_ATTEN > 0) ? atten_of_doca(titer->getD()) : 1.;
            if (q4) {
               // q4: plateau efficiency of the wire x turn-on of its plane vs Qw of this track in this chamber
               int L = 3 * (witer->getModule() - 1) + witer->getLayer(), w = witer->getWire();
               if (L < 1 || L > 24 || w < 1 || w > 96) continue;
               auto it = q4_qsum.find(titer->getItrack());
               double qw = (it == q4_qsum.end()) ? 0. : 0.5 * (it->second.first + it->second.second);
               if (!gDRandom.DecideToAcceptHit(fdc_config->q4_eff[L-1][w-1] * fdc_config->Q4TurnOn(L, qw)))
                  continue;
            }
            else if (titer->getDE() * gain * att * fluct_of_track(titer->getItrack()) * 1e6 < wire_thr_kev)
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

//-----------
// SetRestQ (q4): strip q = FDC:REST_ASCALE x peakAmp for every FDC strip hit with an amplitude, after the random-trigger merge
//-----------
void FDCSmearer::SetRestQ(hddm_s::HDDM *record)
{
   if (fdc_config_t::FDC_REST_ASCALE <= 0) return;
   hddm_s::FdcCathodeHitList hits = record->getFdcCathodeHits();
   for (auto h = hits.begin(); h != hits.end(); ++h)
      if (h->getFdcDigihits().size() > 0)
         h->setQ(fdc_config_t::FDC_REST_ASCALE * h->getFdcDigihit().getPeakAmp());
}

