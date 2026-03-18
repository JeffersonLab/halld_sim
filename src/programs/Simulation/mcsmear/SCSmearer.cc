#include "SCSmearer.h"
#include "DANA/DEvent.h"


//-----------
// sc_config_t  (constructor)
//-----------
sc_config_t::sc_config_t(const std::shared_ptr<const JEvent>& event) 
{
	// default values
	START_SIGMA           = 0.0; // 300ps
	START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200
	START_PADDLE_THRESHOLD  = 0.0;

    // Get the geometry
    DGeometry* locGeometry = DEvent::GetDGeometry(event);

    // Get start counter geometry
    unsigned int MAX_SECTORS=0;
    if (locGeometry->GetStartCounterGeom(sc_pos, sc_norm))  {
		MAX_SECTORS = sc_pos.size();
        for(unsigned int sc_index=0; sc_index<sc_pos.size(); sc_index++) {
            SC_START_Z.push_back( sc_pos[sc_index][0].z() );
        }
            
		double theta = sc_norm[0][sc_norm[0].size() - 2].Theta();
		START_ANGLE_CORR = 1./cos(M_PI_2 - theta);
    }

	// Load data from CCDB
    cout << "Get START_COUNTER/start_parms parameters from CCDB..." << endl;
    map<string, double> startparms;
    if(DEvent::GetCalib(event, "START_COUNTER/start_parms", startparms)) {
		jerr << "Problem loading START_COUNTER/start_parms from CCDB!" << endl;
	} else {
     	START_SIGMA = startparms["START_SIGMA"] ;
     	START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];
	}
	
	cout<<"get START_COUNTER/paddle_mc_efficiency from calibDB"<<endl;
    if(DEvent::GetCalib(event, "START_COUNTER/paddle_mc_efficiency", paddle_efficiencies)) {
    	jerr << "Problem loading START_COUNTER/paddle_mc_efficiency from CCDB!" << endl;
    }

    // Start counter individual paddle resolutions
    vector< vector<double> > sc_paddle_resolution_params;
    if(DEvent::GetCalib(event, "START_COUNTER/TRvsPL", sc_paddle_resolution_params))
        jout << "Error in loading START_COUNTER/TRvsPL !" << endl;
    else {
        if(sc_paddle_resolution_params.size() != MAX_SECTORS)
            jerr << "Start counter paddle resolutions table has wrong number of entries:" << endl
                 << "  loaded = " << sc_paddle_resolution_params.size()
                 << "  expected = " << MAX_SECTORS << endl;

        for(unsigned int i=0; i<MAX_SECTORS; i++) {
            SC_SECTION1_P0.push_back( sc_paddle_resolution_params[i][0] ); 
            SC_SECTION1_P1.push_back( sc_paddle_resolution_params[i][1] );
            SC_BOUNDARY1.push_back( sc_paddle_resolution_params[i][2] );
            SC_SECTION2_P0.push_back( sc_paddle_resolution_params[i][3] ); 
            SC_SECTION2_P1.push_back( sc_paddle_resolution_params[i][4] );
            SC_BOUNDARY2.push_back( sc_paddle_resolution_params[i][5] );
            SC_SECTION3_P0.push_back( sc_paddle_resolution_params[i][6] ); 
            SC_SECTION3_P1.push_back( sc_paddle_resolution_params[i][7] );
            SC_BOUNDARY3.push_back( sc_paddle_resolution_params[i][8] );
            SC_SECTION4_P0.push_back( sc_paddle_resolution_params[i][9] ); 
            SC_SECTION4_P1.push_back( sc_paddle_resolution_params[i][10] );
        }
    }

    map<string,double> sc_mc_correction_factors;
    if(DEvent::GetCalib(event, "START_COUNTER/mc_time_resol_corr", sc_mc_correction_factors)) {
        jout << "Error in loading START_COUNTER/mc_time_resol_corr !" << endl;
    } else {
        SC_MC_CORRECTION_P0 = sc_mc_correction_factors["P0"];
        SC_MC_CORRECTION_P1 = sc_mc_correction_factors["P1"];
    }

}

//------------------------
// GetPaddleTimeResolution
//------------------------
double sc_config_t::GetPaddleTimeResolution(int sector, double sc_local_z)  { 
	double time_resolution = 0.;  // units of ns
	
	// the new 4-region piecewise parameterization is in terms of the pathlength along the paddle
	double dpath = 0.;

	// calculate in local coordinates
	double sc_pos_soss = sc_pos[sector][0].z();   // Start of straight section
	double sc_pos_eoss = sc_pos[sector][1].z() - sc_pos_soss;   // End of straight section
	double sc_pos_eobs = sc_pos[sector][sc_pos[sector].size() - 2].z() - sc_pos_soss;  // End of bend section

	// Calculate hit distance along scintillator relative to upstream end
	if (sc_local_z <= sc_pos_eoss)  // Check to see if hit occured in the straight section
		dpath = sc_local_z ;
	else if(sc_local_z > sc_pos_eoss && sc_local_z <= sc_pos_eobs) //check if in bend section: if so, apply corrections
		dpath = (sc_local_z - sc_pos_eoss)*START_ANGLE_CORR + sc_pos_eoss;
	else // nose section: apply corrections
		dpath = (sc_local_z - sc_pos_eoss)*START_ANGLE_CORR + sc_pos_eoss;

	// look up resolution
	if(dpath < SC_BOUNDARY1[sector]) {
		time_resolution = SC_SECTION1_P0[sector] + SC_SECTION1_P1[sector]*dpath;
	} else if(dpath < SC_BOUNDARY2[sector]) {
		time_resolution = SC_SECTION2_P0[sector] + SC_SECTION2_P1[sector]*dpath;
	} else if(dpath < SC_BOUNDARY3[sector]) {
		time_resolution = SC_SECTION3_P0[sector] + SC_SECTION3_P1[sector]*dpath;
	} else {
		time_resolution = SC_SECTION4_P0[sector] + SC_SECTION4_P1[sector]*dpath;
	}
	
	// If these resolutions come from data, apply correction factors to remove any other contributions
	// note that the factors were determined in units of ps, so we need to convert to and from ns.
	time_resolution = ( (time_resolution*1000. - SC_MC_CORRECTION_P0) / SC_MC_CORRECTION_P1 ) / 1000.;

	// convert ps to ns
	//time_resolution /= 1000.;
	//cout << "z hit = " << sc_local_z << "  dpath = " << dpath << "  time resolution = " << time_resolution << endl;
	return time_resolution;
}


//-----------
// SmearEvent
//-----------
void SCSmearer::SmearEvent(hddm_s::HDDM *record)
{
   hddm_s::StcTruthPointList truthPoints = record->getStcTruthPoints();
        
   hddm_s::StcPaddleList pads = record->getStcPaddles();
   hddm_s::StcPaddleList::iterator iter;
   for (iter = pads.begin(); iter != pads.end(); ++iter) {
      iter->deleteStcHits();
      hddm_s::StcTruthHitList thits = iter->getStcTruthHits();
      hddm_s::StcTruthHitList::iterator titer;
      for (titer = thits.begin(); titer != thits.end(); ++titer) {
      	 // correct simulation efficiencies 
		 if(config->APPLY_EFFICIENCY_CORRECTIONS
		 	&& !gDRandom.DecideToAcceptHit(sc_config->GetEfficiencyCorrectionFactor(iter->getSector())))
		 	continue;

         // smear the time
         hddm_s::StcTruthPointList::iterator piter = FindMatchingTruthPoint(titer, truthPoints);
         // calculate a z-depending timing resolution
         // z is measured from the readout end of the paddles
         double z_pos = 30.;    // default value in the middle, in case we can't find a good point.  this shouldn't happen, but you never know...
         if( piter != truthPoints.end() )
             z_pos = piter->getZ() - sc_config->SC_START_Z[iter->getSector()-1];
    
         double t = titer->getT();
         double NewE = titer->getDE();
         if(config->SMEAR_HITS) {
         	t += gDRandom.SampleGaussian(sc_config->GetPaddleTimeResolution(iter->getSector()-1, z_pos));
         	// smear the energy
         	double npe = titer->getDE() * 1000. *  sc_config->START_PHOTONS_PERMEV;
         	npe = npe +  gDRandom.SampleGaussian(sqrt(npe));
         	NewE = npe/sc_config->START_PHOTONS_PERMEV/1000.;
         }
         if (NewE > sc_config->START_PADDLE_THRESHOLD) {
            hddm_s::StcHitList hits = iter->addStcHits();
            hits().setT(t);
            hits().setDE(NewE);
         }
      }

      if (config->DROP_TRUTH_HITS)
         iter->deleteStcTruthHits();
   }
   if (config->DROP_TRUTH_HITS) {
      hddm_s::StartCntrList stcs = record->getStartCntrs();
      if (stcs.size() > 0)
         stcs().deleteStcTruthPoints();
   }
}

// ----------------------
// FindMatchingTruthPoint
// ----------------------
hddm_s::StcTruthPointList::iterator SCSmearer::FindMatchingTruthPoint(hddm_s::StcTruthHitList::iterator hiter, hddm_s::StcTruthPointList &truthPoints) 
{
    // Match the StcTruthHit with the most likely corresponding StcTruthPoin
    // This is needed since StcTruthHits correspond to detector hits, and so only have time and
    // energy values.   If we want to do something with a z-dependence, e.g. time resolutions,
    // we need the StcTruthPoint, which has a location in detector coordinates.
    // The only thing they have in common in the energy deposited in the scintillator paddles
    // since the StcTruthHit has a propagation time correction applied, so we use that
    // to disambiguate multiple hits in the same paddle
    hddm_s::StcTruthPointList::iterator piter;
    hddm_s::StcTruthPointList::iterator best_piter = truthPoints.end();
    double best_match_deltaE = 100.;
    for( piter = truthPoints.begin(); piter != truthPoints.end(); piter++) {
        if( hiter->getSector() == piter->getSector() ) {
            double deltaE = fabs(hiter->getDE() - piter->getDEdx());
            if(deltaE < best_match_deltaE) {
                best_piter = piter;
                best_match_deltaE = deltaE;
            }
        }
    }

    return best_piter;
}
