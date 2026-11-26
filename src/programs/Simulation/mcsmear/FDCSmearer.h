// Smearing class for forward drift chamber (FDC)

#ifndef _FDCSMEARER_H_
#define _FDCSMEARER_H_

#include "Smearer.h"


class fdc_config_t 
{
  public:
	fdc_config_t(const std::shared_ptr<const JEvent>& event);

	double FDC_TDRIFT_SIGMA;
	double FDC_CATHODE_SIGMA;
	double FDC_PED_NOISE;         // in pC calculated in above
	double FDC_THRESHOLD_FACTOR;  // number of pedestal sigmas for determining sparcification threshold
	//double FDC_HIT_DROP_FRACTION = 0.0; // 1000.0E-9
	double FDC_TIME_WINDOW;
	double FDC_THRESH_KEV;        // fdc anode discriminator threshold
	double FDC_ADC_ASCALE;        // /FDC/digi_scales FDC_ADC_ASCALE; strip peakAmp = q / FDC_ADC_ASCALE, as in data (DFDCHit_factory: q = a_scale * pulse_height)

    static bool FDC_EFFVSDOCA;
    static double FDC_DRIFT_SMEAR_NS;        // FDC:DRIFT_SMEAR_NS, Gaussian drift-time smear (ns) common to the wire and strip hits of one track in one chamber (0 = off)
    static double FDC_CATHODE_SIGMA_SET;     // FDC:CATHODE_SIGMA, override of FDC_CATHODE_SIGMA in microns (< 0 = use CCDB)
    static double FDC_GAIN_SCALE;            // FDC:GAIN_SCALE, gas-gain factor applied to the strip charge q and the wire dE before the thresholds (1 = off)
    static double FDC_GAIN_PKG[4];           // FDC:GAIN_PKG1..4, extra gain factor per FDC package (modules 1-2, 3-4, 5-6, 7-8), times FDC_GAIN_SCALE (1 = off)
    static double FDC_STRIP_THR_SET;         // FDC:STRIP_THRESHOLD, strip charge threshold in pC (< 0 = FDC_THRESHOLD_FACTOR * FDC_PED_NOISE from CCDB)
    static double FDC_GAIN_SIGMA;            // FDC:GAIN_SIGMA s, avalanche gain fluctuation: one log-normal factor exp(s G - s^2/2) per chamber and track, shared by its wire hit and the strips it induces (0 = off)
    static double FDC_DOCA_ATTEN;            // FDC:DOCA_ATTEN a, signal (wire dE and the track's strip q in the chamber) scaled by max(0.05, 1 - a doca^2), doca in cm (0 = off)
    static double FDC_WIRE_THR_KEV_SET;      // FDC:WIRE_THRESHOLD_KEV, anode threshold on the gain-scaled wire dE in keV (< 0 = FDC_THRESH_KEV from CCDB)
    static double FDC_EFFVSDOCA_PAR_SET[4];  // FDC:EFFVSDOCA_PAR0..3, wire efficiency vs doca (defaults = the former hard-coded values)
    static double FDC_RAW_AMP_THR;        // FDC:RAW_AMP_THRESHOLD in ADC counts; < 0 (default): per strip from CCDB /FDC/packageN/strip_thresholds,
                                          // FDC_RAW_AMP_THR_DEFAULT where a table is absent; >= 0: this value for every strip (0: no cut)
    static const double FDC_RAW_AMP_THR_DEFAULT;
    static double FDC_PEAKAMP_ASCALE;     // FDC:PEAKAMP_ASCALE; > 0: strip peakAmp = q / this, the same in every run (the readout threshold
                                          // stays on q / FDC_ADC_ASCALE of the run); <= 0 (default): peakAmp = q / FDC_ADC_ASCALE
    double FDC_EFFVSDOCA_PAR[4];

	// q4 (rho0study-10-2026, 2026-11-20): per-channel FDC response measured from raw hits only (fdccal plugin, BOR record).
	// FDC:Q4_TABLE = path of a text table for the run; empty (default) = off. Lines (# = comment):
	//   S <cp> <strip> <gain> <thr>  cathode plane cp = 2*(gLayer-1) + 1 (plane 1) or + 2 (plane 3), strip 1-192: strip ADC = gain x q (pC),
	//                                readout threshold thr in ADC counts (hit dropped if ADC < thr)
	//   W <L> <wire> <eff>           wire plane L = gLayer 1-24, wire 1-96: efficiency on the charge plateau
	//   T <L> <q50> <width>          wire plane L: turn-on 1/(1 + exp(-(Qw - q50)/width)) vs Qw = (Q_U + Q_V)/2, the mean of the strip ADC
	//                                sums the same track induces on the two cathodes of the chamber (as Qw of the fdccal plugin)
	// With the table: strip ADC (fdcDigihit peakAmp) = gain x q, the readout threshold per strip, the wire hit kept with probability
	// eff x turn-on(Qw); FDC:GAIN_SCALE, GAIN_PKG, RAW_AMP_THRESHOLD, PEAKAMP_ASCALE and the wire keV threshold are not used.
	static std::string FDC_Q4_TABLE;
	// FDC:REST_ASCALE a > 0: after the random-trigger merge, every FDC strip hit gets q = a x peakAmp (MC and random hits alike), so that q
	// is on the scale of the data REST the simulation is compared with (a = /FDC/digi_scales of that REST); <= 0 (default) = off
	static double FDC_REST_ASCALE;
	// q5 (rho0study-10-2026, 2026-11-23): FDC:STRIP_TIME_TABLE = path of a text table of the strip time vs strip amplitude, measured from
	// data (wire t - strip t of isolated cluster pairs, relative to the strips above 800 ADC); empty (default) = off. Lines (# = comment):
	//   <lo> <hi> <centre> <shift> <width>   ADC, ADC, ADC, ns, ns
	// In q4 mode every strip hit gets t += shift(A) + Gaussian(width(A)), A = its ADC, interpolated linearly in ln A between the bin centres
	// and constant beyond the first and last centre, before the readout time window is applied.
	static std::string FDC_STRIP_TIME_TABLE;
	// q5 item 5 (2026-11-23): FDC:STRIP_T0_PKG1..4 (ns, default 0): calibration residue of the strip times of the data REST the simulation is
	// compared with, per package: subtracted from every strip time in q4 mode (data strips come earlier relative to the wire by this much)
	static double FDC_STRIP_T0_PKG[4];
	std::vector<double> stt_lnA, stt_shift, stt_width;
	// quantile form of the table (lines "Q <centre ADC> <q_0 .. q_N>" at p = i/N): the strip time added is drawn from the data's distribution
	// by inverse CDF, the quantiles interpolated linearly in ln A between the bin centres (2026-11-23, to give the data's non-Gaussian tails)
	std::vector<std::vector<double> > stt_q;
	double StripTimeDraw(double adc, double u) const {
		double x = log(adc > 1. ? adc : 1.); size_t n = stt_lnA.size(); size_t i = 0; double f = 0.;
		if (x <= stt_lnA[0]) { i = 0; f = 0.; } else if (x >= stt_lnA[n-1]) { i = n - 2; f = 1.; }
		else { i = 1; while (stt_lnA[i] < x) i++; i--; f = (x - stt_lnA[i]) / (stt_lnA[i+1] - stt_lnA[i]); }
		if (n == 1) { i = 0; f = 0.; }
		auto qd = [&](const std::vector<double> &q) { size_t m = q.size() - 1; double y = u * m; size_t k = (size_t)y; if (k >= m) k = m - 1;
			double g = y - k; return q[k] + g * (q[k+1] - q[k]); };
		double a = qd(stt_q[i]); return (n == 1 || f == 0.) ? a : (1. - f) * a + f * qd(stt_q[i+1]); }
	void StripTime(double adc, double &shift, double &width) const {
		double x = log(adc > 1. ? adc : 1.); size_t n = stt_lnA.size();
		if (x <= stt_lnA[0]) { shift = stt_shift[0]; width = stt_width[0]; return; }
		if (x >= stt_lnA[n-1]) { shift = stt_shift[n-1]; width = stt_width[n-1]; return; }
		size_t i = 1; while (stt_lnA[i] < x) i++;
		double f = (x - stt_lnA[i-1]) / (stt_lnA[i] - stt_lnA[i-1]);
		shift = stt_shift[i-1] + f * (stt_shift[i] - stt_shift[i-1]); width = stt_width[i-1] + f * (stt_width[i] - stt_width[i-1]); }
	bool q4;
	double q4_gain[48][192], q4_thr[48][192], q4_eff[24][96], q4_q50[24], q4_width[24];
	static int q4_cplane(int module, int layer, int plane) { return 2 * (3 * (module - 1) + layer - 1) + ((plane == 3)? 2 : 1); }
	double Q4TurnOn(int L, double qw) const { return 1. / (1. + exp(-(qw - q4_q50[L-1]) / q4_width[L-1])); }

	vector< vector<double > > channel_efficiencies;
	vector< vector<double > > channel_raw_amp_thresholds;   // same row order as channel_efficiencies; empty rows for the anode wires
	
	double GetEfficiencyCorrectionFactor(hddm_s::FdcCathodeStripList::iterator &siter) {
		// cathode strips
		int gPlane = 9*(siter->getModule()-1) + 3*(siter->getLayer()-1) 
             		+ (siter->getPlane()-1);               // starts counting at 0
        int element = siter->getStrip();

        //cout << "module = " << siter->getModule() << " layer = " << siter->getLayer()
        //   << " plane = " << siter->getPlane() << " gPlane = " << gPlane 
        //   << " strip = " << siter->getStrip() << endl;

		return channel_efficiencies.at(gPlane).at(element-1);
	}

	double GetRawAmpThreshold(hddm_s::FdcCathodeStripList::iterator &siter) {
		// readout threshold of a cathode strip on raw_amp = q / FDC_ADC_ASCALE, in ADC counts
		if (FDC_RAW_AMP_THR >= 0)
			return FDC_RAW_AMP_THR;
		int gPlane = 9*(siter->getModule()-1) + 3*(siter->getLayer()-1) 
             		+ (siter->getPlane()-1);               // starts counting at 0
		int element = siter->getStrip();
		const vector<double> &row = channel_raw_amp_thresholds.at(gPlane);
		return (element >= 1 && element <= (int)row.size())? row[element-1] : FDC_RAW_AMP_THR_DEFAULT;
	}

	double GetEfficiencyCorrectionFactor(hddm_s::FdcAnodeWireList::iterator &witer) {
		// anode wires
		int gPlane = 9*(witer->getModule()-1) + 3*(witer->getLayer()-1) + 1;  // starts counting at 0
        int element = witer->getWire();

        //cout << "module = " << witer->getModule() << " layer = " << witer->getLayer()
        //  //<< " plane = " << witer->getPlane() 
        //   << " gPlane = " << gPlane 
        //   << " wire = " << witer->getWire() << endl;

		return channel_efficiencies.at(gPlane).at(element-1);
    }

	double GetEfficiencyVsDOCA(double doca) {
        // Introduce a DOCA-dependent efficiency that cuts off the hit efficiency
        // for tracks that cross through the cell at a large distance from the wire.
        double eff = FDC_EFFVSDOCA_PAR[0] - FDC_EFFVSDOCA_PAR[1] *
                     (1 / (pow(FDC_EFFVSDOCA_PAR[2] - doca, 2) +
                           pow(FDC_EFFVSDOCA_PAR[3], 2))
                      - 2 * FDC_EFFVSDOCA_PAR[2] * doca /
                          pow(pow(FDC_EFFVSDOCA_PAR[2], 2) +
                              pow(FDC_EFFVSDOCA_PAR[3], 2), 2)
                     );
        return (!FDC_EFFVSDOCA)? 1 : (eff < 0)? 0 : eff*eff;
	}
};



class FDCSmearer : public Smearer
{
  public:
	FDCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		fdc_config = new fdc_config_t(event);
	}
	~FDCSmearer() {
		delete fdc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);
	static void SetRestQ(hddm_s::HDDM *record);   // FDC:REST_ASCALE, called after the random-trigger merge

  private:
  	fdc_config_t  *fdc_config;
};



#endif // _FDCSMEARER_H_
