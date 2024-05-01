// Smearing class for start counter (SC)

#ifndef _SCSMEARER_H_
#define _SCSMEARER_H_

#include "Smearer.h"
#include "TMath.h"


class sc_config_t 
{
  public:
	sc_config_t(const std::shared_ptr<const JEvent>& event);

    double GetPaddleTimeResolution(int sector, double sc_local_z);

	double GetEfficiencyCorrectionFactor(int sector) {
		return paddle_efficiencies.at(sector-1);
	}
	
	double START_SIGMA;
	double START_PHOTONS_PERMEV;
	double START_PADDLE_THRESHOLD;
	
	double START_ANGLE_CORR;

	vector<double> paddle_efficiencies;

    // Start counter geometry parameters
    vector<vector<DVector3> >sc_norm; 
    vector<vector<DVector3> >sc_pos;
    vector<double> SC_START_Z;

    // Start counter resolution parameters
    vector<double> SC_BOUNDARY1, SC_BOUNDARY2, SC_BOUNDARY3;
    vector<double> SC_SECTION1_P0, SC_SECTION1_P1;
    vector<double> SC_SECTION2_P0, SC_SECTION2_P1;
    vector<double> SC_SECTION3_P0, SC_SECTION3_P1;
    vector<double> SC_SECTION4_P0, SC_SECTION4_P1;

    double SC_MC_CORRECTION_P0, SC_MC_CORRECTION_P1;

};


class SCSmearer : public Smearer
{
  public:
	SCSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
		sc_config = new sc_config_t(event);
	}
	~SCSmearer() {
		delete sc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

  protected:
    hddm_s::StcTruthPointList::iterator FindMatchingTruthPoint(hddm_s::StcTruthHitList::iterator hiter, hddm_s::StcTruthPointList &truthPoints);
	
  private:
  	sc_config_t  *sc_config;
};



#endif // _SCSMEARER_H_
