// Classes to store configuration information for mcsmear

#ifndef _MCSMEAR_CONFIG_H_
#define _MCSMEAR_CONFIG_H_

#include "units.h"
#include <DANA/DEvent.h>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <DVector3.h>

using std::string;
using std::vector;
using std::cout;

// #include <DANA/DApplication.h>
#include "DRandom2.h"
#include <JANA/JEvent.h>



// external function definitions from SampleGaussian.cc
double SampleGaussian(double sigma);
double SamplePoisson(double lambda);
double SampleRange(double x1, double x2);


// Overall configuration parameters
class mcsmear_config_t 
{
  public:
	mcsmear_config_t();
	~mcsmear_config_t();

	//-----------
	// SetSeeds
	//-----------
	void SetSeeds(const char *vals);

	// member variables
	bool ADD_NOISE;
	bool DROP_TRUTH_HITS;
	bool SMEAR_HITS;
	bool DUMP_RCDB_CONFIG;
	bool SKIP_READING_RCDB;
	bool MERGE_TAGGER_HITS;

	//bool SMEAR_BCAL;
	//bool FDC_ELOSS_OFF;
	bool IGNORE_SEEDS;
	double TRIGGER_LOOKBACK_TIME;
	bool APPLY_EFFICIENCY_CORRECTIONS;
	bool APPLY_HITS_TRUNCATION;

    bool FCAL_ADD_LIGHTGUIDE_HITS;
    double FCAL_LIGHTGUIDE_SCALE_FACTOR;
	
	bool FCAL_NEW_TIME_SMEAR;
	// flags to pass command line info to subdetector classes
	double BCAL_NO_T_SMEAR;
	double BCAL_NO_DARK_PULSES;
	double BCAL_NO_SAMPLING_FLUCTUATIONS;
	double BCAL_NO_SAMPLING_FLOOR_TERM;
	double BCAL_NO_POISSON_STATISTICS;
	double BCAL_NO_FADC_SATURATION;
	double BCAL_NO_SIPM_SATURATION;

	// list of detectors with hits to smear
	std::string DETECTORS_TO_LOAD="all";
	
	
#ifdef HAVE_RCDB
    void LoadRCDBConnection();
	bool ParseRCDBConfigFile(int runNumber);
#endif  // HAVE_RCDB

    std::map<std::string, std::map<std::string, double> > readout;
};


#endif  // _MCSMEAR_CONFIG_H_
