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

    static double FDC_DRIFT_SMEAR_NS;        // FDC:DRIFT_SMEAR_NS, Gaussian drift-time smear (ns) common to the wire and strip hits of one track in one chamber (0 = off)
    static double FDC_CATHODE_SIGMA_SET;     // FDC:CATHODE_SIGMA, override of FDC_CATHODE_SIGMA in microns (< 0 = use CCDB)
    static double FDC_GAIN_SCALE;            // FDC:GAIN_SCALE, gas-gain factor applied to the strip charge q and the wire dE before the thresholds (1 = off)
    static double FDC_GAIN_PKG[4];           // FDC:GAIN_PKG1..4, extra gain factor per FDC package (modules 1-2, 3-4, 5-6, 7-8), times FDC_GAIN_SCALE (1 = off)
    static double FDC_STRIP_THR_SET;         // FDC:STRIP_THRESHOLD, strip charge threshold in pC (< 0 = FDC_THRESHOLD_FACTOR * FDC_PED_NOISE from CCDB)
    static double FDC_GAIN_SIGMA;            // FDC:GAIN_SIGMA s, avalanche gain fluctuation: one log-normal factor exp(s G - s^2/2) per chamber and track, shared by its wire hit and the strips it induces (0 = off)
    static double FDC_DOCA_ATTEN;            // FDC:DOCA_ATTEN a, signal (wire dE and the track's strip q in the chamber) scaled by max(0.05, 1 - a doca^2), doca in cm (0 = off)
    static double FDC_WIRE_THR_KEV_SET;      // FDC:WIRE_THRESHOLD_KEV, anode threshold on the gain-scaled wire dE in keV (< 0 = FDC_THRESH_KEV from CCDB)

	vector< vector<double > > channel_efficiencies;
	
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
	
  private:
  	fdc_config_t  *fdc_config;
};



#endif // _FDCSMEARER_H_
