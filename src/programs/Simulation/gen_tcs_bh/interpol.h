#define interpol_h
#include "TabUtils.h"
#include "TableProvider.h"
#include "TreatOptions.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

double linear_interpol_tcs5D( 
	TableProvider ddvcs_table_provider,float &cross_BH,float &cross_TCS,
        float &polpp, float &polpm, float &polmp, float &polmm,
        double xbj, double mt, double Q2, double PhiCM, double ThetaCM,
        float xbj_grid[], float mt_grid[],float Q2_grid[], float ThetaCM_grid[], float PhiCM_grid[],
        int i, int j, int k, int n, int o);

double linear_interpol_tcs6D_4( 
	TableProvider ddvcs_table_provider,
	float &W_BH, float &W_TCS, float &pp,
	float &pm , float &mp , float &mm , float &bha, float &bhb,  
	double mt, double Q2, double PhiCM, double ThetaCM,
	float mt_grid[],float Q2_grid[], float ThetaCM_grid[], float PhiCM_grid[],
	int targetpoldir, int beampoltype,
	int i, int j, int k, int n, int o,int l
);


double getval_tcs(Vars_Sef_tcs sef,int crA);
double getval_tcs_pol(Vars_Sef_tcs_pol sef,int crA);


// more options unused, see in ~/Generator/Archives/Generator_grid_Sept2018_save/includes/interpol.h
// extrapol and other interpol options
