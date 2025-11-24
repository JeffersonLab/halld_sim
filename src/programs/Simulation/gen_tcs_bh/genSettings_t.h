#ifndef genSettings_t_h
#define genSettings_t_h
#include "UTILITIES/MyReadConfig.h"
#include <iostream>
#include <string>
using namespace std;

// CONFIGURATION SETTINGS FOR TCS GENERATOR
struct genSettings_t {
 	int reaction;       //reaction (1 -> TCS; 11-> TCS phase space)
						//    1) tcs = TCS exclusive production of lepton pair from gamma or e beam"<<endl;
					   	//	  - phase space for exclusive processes (flags):"<<endl;
					    //	  11) ps_eephoto_fix = TCS phase space for fix target exp. (no weight),"<<endl;

    // --- Beam configuration -------------------------------------------------
    int beamtype;          // Electron beam=1 (equivalent photon approximation) or photon beam=0

    double EphotonMin;     // Energy range of the photon (maximal range: 5 < Egamma < 11.5 GeV or less than electron energy)
    double EphotonMax;     // photon energy maximum
    double Eelectron;      // In case of initial electron beam, energy of the incoming electron
						   // In case of a fix energy photon beam, electron energy before radiator (used to calculate the polarization)


    double thetaphotoMax;   // In case of initial electron beam, or collimated real photon (not only along z-axis), polar angle limit of the photon (real and/or quasi-real). -> recommended to let it at 0
							// Default = 0 for real photon. If set at 0, the beam is along axis without any dispersion.
							// If set at 0 for quasi-real beam, no cut is applied. Average angles without cut: 5e-5 rad for bremsstrahlung, 1e-2 rad for quasi-real.

    // --- Event generation ----------------------------------------------------
    int nTotEvents;        	// number of events to generate
    int runNum;             //run number
    int outLepton;        	// Outgoing particle to study: 1=electron, 2=muon
							// Warning: mass correction is included only for kinematic variables, not in cross sections

    // --- Target configuration ------------------------------------------------
    double targetLength;  	// target length (cm); Target lenght in cm and target atoms {A,N} or choose an element instead of A (for LH2 target, set (1,1) or (1001,1)) -> relevant only if use of quasi-real photon beam
							// target type and code = LH2: 1001; He: 1002; LiH: 1003; Carbon: 1006; NH3: 1007;  dry air: 1011; Al: 1013; Fe: 1016; Pb: 1082; lead-glas: 1083.
							//	or with atoms, use numbers H (1), He (2), Be (4), N-H3 (7) , Au (79).

    int    A_target;      // target atomic
    int    Z_target;      // target proton

    int protonOrNeutron;  // TCS on Proton (1) or TCS on Neutron (2)

    // --- Polarization --------------------------------------------------------
    double polBeamDeg;    	// Beam polarization: degree of linear e- or circular gamma polarization in lab frame (<1)
							// If electron beam, will calculate P_circ as a function of the photon energy
							// If circularly polarized photon beam from bremstrahlung: calculated as a function of energy
							// If linearly polarized photon: set average polarization (energy dependence of P is not included)

    int    beamPolarType; 	// Real photon polarization type: set 0 if no polarizer used (circular polarization rate is calculated as a fonction of energy transfered to the photon, or it is fixed to beam polarization in case of fix photon energy)
							// set 0 if a linearly polarized electron beam is used or a circularly polarized photon beam (default option)
							// if a linear polarizer is used for photons: 1 = along x-axis, 2=along y-axis, 3=45^0

    int    targetPolDir;  	// Target nucleon polarization, incoherent scattering on proton
							// 0=unpolarized, +1 or +2 = transverse (perp to beam: 1 = x-axis, 2 = y-axis), +3 = longitudinal (same orientation as beam)

    double polTargetDeg;    // degree of target polarization in lab frame(between -1 and +1)

    // --- Kinematic limits ----------------------------------------------------
    double mt_Min;        // minimal -t; limits in -t (maximal range: 0.04 < -t < 2.02 GeV2)
    double mt_Max;        // maximal -t

    double Qp2Min;        // minimal Q'²; limits in Q'2 (maximal range: 3.8 < Q'2 < 9.2 GeV2)
    double Qp2Max;        // maximal Q'²

    double thetaCMMin;    // outgoing lepton CM theta_min (deg); outgoing lepton angle in CM frame (maximal range= 30 < th < 150 degree)
    double thetaCMMax;    // outgoing lepton CM theta_max (deg)

    double Q2Max;         // maximal allowed Q² in case of electron beam, maximal Q2 allowed (max=0.3 GeV2)

    // --- Radiative corrections -----------------------------------------------
    int radCorrType;      //LO (0) or include radiative corrections:
							//(1) only external radiation of beam (2) + internal real corrections (3) + virtual corrections


    double eCut;          // if radiative corrections apply, set energy cut-off (max) . applied on electron beam and leptons out
							// Warning: Ecut has to be less than Eelectron - 1.5 GeV

    // --- Output options ------------------------------------------------------
    int outFormat;        // HDDM is default
							// Print in ROOT only: 0, ROOT+HEP: 1, HEP only: 2, input for SIMC only: 3, SIMC+ROOT:4
							// simc text only, angular and momentum cuts (below): 5
							// simc+root, angular and momentum cuts (below): 6

    // --- SIMC-specific parameters --------------------------------------------
    double thetaHMS;      // If for SIMC, theta HMS vs beam axis (deg) only with option SIMC (opt 5 or 6) +-25%
							// default=0 if not used
    double thetaSHMS;     // If for SIMC, theta SHMS vs beam axis (deg) only with option SIMC (opt 5 or 6) +(25%)
							// default=0 if not used

    double pHMS;          // If for SIMC, momentum HMS (opt 5 or 6) +-25%
							// default=0 if not used
    double pSHMS;         // If for SIMC, momentum SHMS (opt 5 or 6) +-25%
							// default=0 if not used

    // --- File configuration --------------------------------------------------
    char outFile[180];    // output filename
    //string outFile;
	int indexrun;
	long long int seedentry;		 // seed for random number generator


// LOAD SETTINGS FROM CONFIG FILE
        genSettings_t() {
        MyReadConfig * ReadFile = new MyReadConfig();
        ReadFile->ReadConfigFile("tcs_bh.cfg");  // generator config file
		reaction    	 = *(ReadFile->GetConfig1Par("reaction"));       // TCS
        // Beam settings
        beamtype        = *(ReadFile->GetConfig1Par("beamtype"));       // fixed-energy photon beam
        EphotonMin      = *(ReadFile->GetConfig1Par("EphotonMin"));
        EphotonMax      = *(ReadFile->GetConfig1Par("EphotonMax"));
        Eelectron       = *(ReadFile->GetConfig1Par("Eelectron"));

        thetaphotoMax   = *(ReadFile->GetConfig1Par("thetaphotoMax"));
        nTotEvents      = *(ReadFile->GetConfig1Par("nTotEvents")); //50;
		runNum          = *(ReadFile->GetConfig1Par("runNum"));
        outLepton       = *(ReadFile->GetConfig1Par("outLepton"));       // electron

        targetLength    = *(ReadFile->GetConfig1Par("targetLength"));    // cm
        A_target        = *(ReadFile->GetConfig1Par("A_target"));    // LH2
        Z_target        = *(ReadFile->GetConfig1Par("Z_target"));
        protonOrNeutron = *(ReadFile->GetConfig1Par("protonOrNeutron"));      // proton

        polBeamDeg      = *(ReadFile->GetConfig1Par("polBeamDeg"));
        beamPolarType   = *(ReadFile->GetConfig1Par("beamPolarType"));       // x-axis
        targetPolDir    = *(ReadFile->GetConfig1Par("targetPolDir"));       // unpolarized
        polTargetDeg    = *(ReadFile->GetConfig1Par("polTargetDeg"));

        mt_Min          = *(ReadFile->GetConfig1Par("mt_Min"));
        mt_Max          = *(ReadFile->GetConfig1Par("mt_Max"));
        Qp2Min          = *(ReadFile->GetConfig1Par("Qp2Min"));
        Qp2Max          = *(ReadFile->GetConfig1Par("Qp2Max"));

        thetaCMMin      = *(ReadFile->GetConfig1Par("thetaCMMin"));
        thetaCMMax      = *(ReadFile->GetConfig1Par("thetaCMMax"));
        Q2Max           = *(ReadFile->GetConfig1Par("Q2Max"));
        radCorrType     = *(ReadFile->GetConfig1Par("radCorrType"));      // LO only
        eCut            = *(ReadFile->GetConfig1Par("eCut"));

        outFormat       = *(ReadFile->GetConfig1Par("outFormat"));
        thetaHMS        = *(ReadFile->GetConfig1Par("thetaHMS"));
        thetaSHMS       = *(ReadFile->GetConfig1Par("thetaSHMS"));
        pHMS            = *(ReadFile->GetConfig1Par("pHMS"));
        pSHMS           = *(ReadFile->GetConfig1Par("pSHMS"));
        // Filenames
		indexrun = *(ReadFile->GetConfig1Par("indexrun"));
		seedentry = *(ReadFile->GetConfig1Par("seedentry")); 	 // default seed
        TString outFile0 = ReadFile->GetConfigName("outFile");
        //strcpy(outFile, ReadFile->GetConfigName("outFile").Data());
        strcpy(outFile, outFile0.Data());
        //strcpy(outFile, "tcs_bh_output.hddm");
    }
};


#endif