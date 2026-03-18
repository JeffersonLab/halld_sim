// Smearing class for forward calorimeter (FCAL)

#ifndef _FCALSMEARER_H_
#define _FCALSMEARER_H_

#include "Smearer.h"

#include <FCAL/DFCALGeometry.h>


class fcal_config_t 
{
  public:
	fcal_config_t(const std::shared_ptr<const JEvent>& event, const DFCALGeometry *fcalGeom, mcsmear_config_t *in_config);

	double FCAL_PHOT_STAT_COEF;
	double FCAL_BLOCK_THRESHOLD;
	double FCAL_TSIGMA;
	vector<double> FCAL_PEDS; 
	vector<double> FCAL_GAINS;
	double FCAL_MC_ESCALE;
        double FCAL_ADC_ASCALE;
	double FCAL_PED_RMS;
        double FCAL_THRESHOLD_SCALING;
        double FCAL_THRESHOLD;
        double FCAL_INTEGRAL_PEAK;
	double FCAL_ENERGY_WIDTH_FLOOR;	
	double INSERT_PHOT_STAT_COEF;
	double INSERT_ENERGY_WIDTH_FLOOR;	
	double FCAL_ENERGY_RANGE;
	bool FCAL_ADD_LIGHTGUIDE_HITS;
    	double FCAL_LIGHTGUIDE_SCALE_FACTOR;
	bool FCAL_NEW_TIME_SMEAR;
 	double FCAL_TIME_A;
  	double FCAL_TIME_B;
  	double FCAL_TIME_C;
  	double FCAL_TIME_D;
  	double FCAL_TIME_E;
	
	vector< vector<double > > block_efficiencies;
	
	double GetEfficiencyCorrectionFactor(double row, double column) {
		return block_efficiencies.at(row).at(column);
	}
};



class FCALSmearer : public Smearer
{
 public:
 FCALSmearer(const std::shared_ptr<const JEvent>& event, mcsmear_config_t *in_config) : Smearer(event, in_config) {
    // Get the FCAL geometry
    event->GetSingle(fcalGeom);

    fcal_config = new fcal_config_t(event, fcalGeom, in_config);
    fcal_config->FCAL_NEW_TIME_SMEAR = in_config->FCAL_NEW_TIME_SMEAR;
    
  }
  ~FCALSmearer() {
    delete fcal_config;
  }
  
  void SmearEvent(hddm_s::HDDM *record);
  
 private:
  fcal_config_t  *fcal_config;
  const DFCALGeometry *fcalGeom;
};


#endif // _FCALSMEARER_H_
