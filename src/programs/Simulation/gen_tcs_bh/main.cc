#include <stdlib.h>
#include "Constants.h"
#include "FormFactors.h"
#include "Options_tcs.h"
#include "TreatOptions.h"
#include "Utils.h"
#include "TabUtils.h"
#include "interpol.h"
#include "TableProvider.h"
#include "PDFparam.h"
#include "TreeInit.h"
#include "Polar.h"
#include "Kinematics.h"
#include "LeptonPair.h"
#include "PartUtils.h"
#include "ReactionKinTwo.h"
#include <chrono>
#include "HddmOut.h" // include after test
#include "genSettings_t.h"


using namespace std;
// https://hallaweb.jlab.org/wiki/index.php/DDVCS_and_TCS_event_generator

int reaction;
int main(int argc, char *argv[]){
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();
	genSettings_t genSettings;
	reaction = genSettings.reaction;
	static int indexrun = genSettings.indexrun;
	int runNum = genSettings.runNum;

	cout << "Program starts" << endl;
    const bool DEBUG_TABLE = true;
    if(DEBUG_TABLE) cout << "[DT] Initialization... " << endl << "[DT] Table Provider creation ... ";
    TableProvider ddvcs_table_provider;

	FillTables filltable;
	Polar pol;

	RadProcess rad_ob;
	LeptonPair pair;
	PartUtils part_op;
	ReactionKinTwo CMlabT;

	long long int seedentry = genSettings.seedentry;
	if (argc > 3) seedentry = (long long int) stoll(argv[3]);
    cout<<"Generate events for reaction: "<<argv[1]<<" (index "<<reaction<<")"<<endl;

	TFile *file;  // Declare the file pointer
	file = new TFile(Form("tcs_bh_output_%d.root", indexrun),"RECREATE");
	if ( file->IsOpen() ) printf("File opened successfully\n");

/*	TFile *file = nullptr;
	if (genSettings.outFormat != 00) {
		file = new TFile(Form("tcs_bh_output_%d.root", indexrun), "RECREATE");
	}
*/
	TTree *SIM_Tree = new TTree("SIM_Tree","generated events");
	TTree *Dump_Tree = new TTree("Dump_Tree","All 4-vectors for 50 first events / initial parameters");

	// HDDM output file
	HddmOut hddmGo(genSettings.outFile);
	/////////////////////////////////////////////////////////////////////////////
//	if (genSettings.outFormat==00){ }else{
	// save in the tree
	static long long int FirstEvent, EventNumber, TrueEventNumber;
	static float param_initfile[30]; // all parameters from input file
	static int rad_event, VirtualFlag, FlagSing, target_spindir, beam_spindir;
	static float epsilon, phase_space, flux_rc;

	static float W_BH, W_TCS, W_tot_unpol, cross_tot_unpol,cross_BH,cross_TCS, BSA,  TSA, BTSA, cross_tot_pol, W_tot_pol, cross_tot_pol_beam, W_tot_pol_beam, cross_tot_pol_target, W_tot_pol_target, cross_BH_x, cross_BH_y, thetamin_nocut; 
	static double Flux_qr, Flux_bmr, theta_beam, ttmin, Psi_s,WW, phi_s, theta_s, CosThetagg, CosThetaggCM; 
	static double ALV_minus_lab[4], ALV_plus_lab[4], ALV_gamma_in[4], ALV_el_in[4], ALV_Recoil_lab[4], ALV_el_out[4];

	// only in the code
	float bremgam_integral, bremgam_integral_internal, targetlenght_rc, d_internal;
	float polpp, polpm, polmm, polmp, phisbin, psisbin, Ebin, mtbin,  Qp2bin,  thetabin, phibin, polarbeaminput; 
	TLorentzVector LV_minus_lab, LV_plus_lab, LV_gamma_lab,LV_Recoil_lab,LV_el_in, LV_el_out, LV_gamma_out_lab, LV_minus_CMeP, LV_plus_CMeP, LV_plus,LV_minus, LV_minus_CMV, LV_plus_CMV, LV_gamma_CMeP, LV_Recoil_CMeP, LV_Recoil, LV_target_lab;
	double theta_gamma, Gamma, Beta, PVirtual_CMeP, EVirtual_CMeP, test, r1, r2, r3, r4;

	int check=0, i1, i2, i3;
    float Egamma_bintable[NEB+1],tt_bintable[NT+1],Qp2_bintable[NQp2+1],phi_bintable[NPhi+1],theta_bintable[NTh+1];
	float  phis_bintable[NPhis+1], psis_bintable[NPsis+1];


	float thmEb[thmE+1]={1000}, thmQb[thmQ+1]={1000}, thmTb[thmT+1]={1000};
    double thm[thmE][thmT][thmQ];

	ifstream cutfile;
	string cfi;
	ofstream hepfile;//(Form("tcs_bh_output_hep_%d.dat", indexrun), std::ofstream::out);
	ofstream kinfile;//(Form("tcs_bh_output_kin_%d.dat", indexrun), std::ofstream::out);
	ofstream logfile(Form("tcs_bh_output_%d.log", indexrun), std::ofstream::out);

	/////////////////////////////////////////////////////////////////////////
	cout<<"runindex: "<<indexrun<<" seed: "<<seedentry<<" time "<<time(NULL)<<endl;

	if (seedentry==0) srand((unsigned int) time(NULL));
	else {
		seedentry/=2;
		seedentry += (unsigned int) time(NULL)/2;
		seedentry += indexrun;
		srand(seedentry);
		}

	/////////////////////////////////////////////////////////////////////////
	printf("\n Init program:  ");


	if (reaction==1 || reaction==11 ){ // tcs

		// binning theta_phi cut
        for (int i=0; i<=thmE; i++){
            thmEb[i]=(float) LinearBins(thmE, thmEmin,thmEmax,i);
			}
		for (int i=0; i<=thmQ; i++){
			thmQb[i]=(float) LinearBins(thmQ, thmQmin,thmQmax,i);
			}
		for (int i=0;i<=thmT;i++){
            thmTb[i]=(float) LinearBins(thmT, thmTmin,thmTmax,i);
			}


		// JLab TCS mode
		if (reaction==1){
			check=Variables_TCS(); // read external option file
			if (check!=1) { cout<<"ERROR: CHECK OPTION FILE AND DON'T CHANGE LINES INSIDE\nRETURN"<<endl; file->Close(); return 0; }
			check=VariablesValues_TCS();  // read options and check validity
			if (check!=1) { cout<<"ERROR: CHECK VARIABLES VALUES IN OPTION FILE\n RETURN"<<endl;file->Close();return 0; }

			printf("\n Load tables... this step may take some time...");
				check = filltable.FillTable(reaction ,ddvcs_table_provider); // read cross sections from external file and put it in arrays
				if (check!=1) { cout<<"ERROR: no tables file \n RETURN"<<endl;file->Close(); return 0; }
				else printf ("..done\n ");

			bremgam_integral = rad_ob.InitBrmProfile(Eelectron,EMINI, EMAXI);

				if (verbose){
					cout<<"Bin E:   "<<EMINI<<" "<<EMAXI<<" " <<NEB<<endl;
					cout<<"Bin -t   "<<TMINI<<" "<<TMAXI<<" "<<NT<<endl;
					cout<<"Bin Qp2: "<<QP2MINI<<" "<<QP2MAXI<<" "<<NQp2<<endl;
					cout<<"Bin phi: "<<PHIMINI<<" "<<PHIMAXI<<" "<<NPhi<<endl;
					cout<<"Bin th:  "<<THMINI<<" "<<THMAXI<<" "<<NTh<<endl;
					}

			// binning cross sections
			for (int i=0; i<=NEB; i++){
				Egamma_bintable[i] = (float) LinearBins(NEB,EMINI,EMAXI,i);
				if (verbose) cout<<" Eb["<<i << "] = "<<Egamma_bintable[i];
				}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NT; i++){
					tt_bintable[i] = (float) LinearBins(NT,TMINI,TMAXI,i);
					if (verbose) cout<<" -t["<<i<<"] = "<<tt_bintable[i];
					}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NQp2; i++){
					Qp2_bintable[i] = (float) LinearBins(NQp2,QP2MINI,QP2MAXI,i);
					if (verbose) cout<<" Qp2["<<i<<"] = "<<Qp2_bintable[i];
					}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NPhi; i++){
					phi_bintable[i] = (float) LinearBins(NPhi,PHIMINI*PI/180.,PHIMAXI*PI/180.,i);
					if (verbose) cout<<" Phi["<<i<<"] = "<<phi_bintable[i];
					}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NTh; i++){
					theta_bintable[i] = (float) LinearBins(NTh,THMINI*PI/180.,THMAXI*PI/180.,i);
					if (verbose) cout<<" Th["<<i<<"] = "<<theta_bintable[i];
					}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NPsis; i++){
					psis_bintable[i] = (float) LinearBins(NPsis,PSISMINI*PI/180.,PSISMAXI*PI/180.,i);
					if (verbose) cout<<" Psis["<<i<<"] = "<<psis_bintable[i];
					}

			if (verbose) cout<<" "<<endl;
				for (int i=0; i<=NPhis; i++){
					phis_bintable[i] = (float) LinearBins(NPhis,PHISMINI*PI/180.,PHISMAXI*PI/180.,i);
					if (verbose) cout<<" phis["<<i<<"] = "<<phis_bintable[i];
					}

			if (verbose) cout<<" "<<endl;

	    // TCS for generic fix target experiment, phase space only
		} else if (reaction==11){

		check = Variables_PSEEPHOTO_FIX();
		if (check!=1) { cout<<"ERROR: CHECK OPTION FILE AND DON'T CHANGE LINES INSIDE  "<<endl; file->Close(); return 0; }
		check=VariablesValues_PSTCSFIX();
		if (check!=1) { cout<<"ERROR: CHECK VARIABLES VALUES IN OPTION FILE   \n RETURN"<<endl; file->Close(); return 0; }
			}

		// BH peaks, read table to position
        cfi=Form("Data/scanBHsing_fullrange.dat");
        cutfile.open(cfi.c_str()); cout<<"table "<<cfi<<endl;
        if (!cutfile) { cout<<"  ERROR: Table for thetaphi cut cannot open. use set.csh "<<endl; file->Close(); return 7;}

		i1=0;i2=0;i3=0;
            while (cutfile.good()){
                if (!(cutfile>>r1)) break;
                    cutfile>>r2>>r3>>r4;
                    thm[i1][i2][i3]=r4;
                    i3++;
                        if (i3==thmQ){
                            i3=0;
                            i2++;
                                if (i2==thmT){
                                    i2=0; i1++;
									}
                                if (i1==thmE) break;
							}
                }
		cutfile.close();

		i1=-1; i2=-1; i3=-1;
		phase_space = (mt_max-mt_min)*(Qp2max-Qp2min) * (theta_max-theta_min)*PI/180.*2.*PI;

		// Tree for TCS events
		SIM_Tree->Branch("ALV_minus_lab",&ALV_minus_lab,"ALV_minus_lab[4]/D");
		SIM_Tree->Branch("ALV_plus_lab",&ALV_plus_lab,"ALV_plus_lab[4]/D");
		SIM_Tree->Branch("ALV_gamma_in",&ALV_gamma_in,"ALV_gamma_in[4]/D");
		SIM_Tree->Branch("ALV_Recoil_lab",&ALV_Recoil_lab,"ALV_Recoil_lab[4]/D");

		SIM_Tree->Branch("Egamma",&Egamma,"Egamma/D");
		SIM_Tree->Branch("Qp2",&Qp2,"Qp2/D");
		SIM_Tree->Branch("tt",&tt,"tt/D");
		SIM_Tree->Branch("ttmin",&ttmin,"ttmin/D");

		SIM_Tree->Branch("Phi_CMV",&Phi_CMV,"Phi_CMV/D");
		SIM_Tree->Branch("Theta_CMV",&Theta_CMV,"Theta_CMV/D");
		if (beamtype!=0) {
			SIM_Tree->Branch("ALV_el_in",&ALV_el_in,"ALV_el_in[4]/D");
			SIM_Tree->Branch("ALV_el_out",&ALV_el_out,"ALV_el_out[4]/D");
			SIM_Tree->Branch("yy",&yy ,"yy/D");
			SIM_Tree->Branch("WW",&WW ,"WW/D");
			SIM_Tree->Branch("Q2",&Q2,"Q2/D");
			SIM_Tree->Branch("epsilon",&epsilon,"epsilon/F");
			SIM_Tree->Branch("VirtualFlag",&VirtualFlag,"VirtualFlag/I");
			SIM_Tree->Branch("Flux_qr",&Flux_qr,"Flux_qr/D");
			SIM_Tree->Branch("Flux_bmr",&Flux_bmr ,"Flux_bmr/D");
			SIM_Tree->Branch("theta_gamma",&theta_gamma ,"theta_gamma/D");
			SIM_Tree->Branch("theta_beam",&theta_beam ,"theta_beam/D");
			}

		if (targetpoldir==1 || targetpoldir==2){
			SIM_Tree->Branch("phi_s", &phi_s, "phi_s/D"); // transverse target polar only
			SIM_Tree->Branch("theta_s", &theta_s, "theta_s/D"); // transverse target polar only
			}
		SIM_Tree->Branch("CosThetagg",&CosThetagg,"CosThetagg/D");
		SIM_Tree->Branch("phi_beam",&phi_beam ,"phi_beam/D");
		if (beampolartype==1) SIM_Tree->Branch("Psi_s", &Psi_s, "Psi_s/D");

		if (reaction<=10){
			if (beamtype!=0){
				SIM_Tree->Branch("cross_tot_unpol",&cross_tot_unpol ,"cross_tot_unpol/F");
				SIM_Tree->Branch("cross_tot_pol", &cross_tot_pol, "cross_tot_pol/F");
				SIM_Tree->Branch("cross_tot_pol_beam", &cross_tot_pol_beam, "cross_tot_pol_beam/F");
				SIM_Tree->Branch("cross_tot_pol_target", &cross_tot_pol_target, "cross_tot_pol_target/F");
				SIM_Tree->Branch("cross_TCS",&cross_TCS ,"cross_TCS/F");
				SIM_Tree->Branch("cross_BH",&cross_BH ,"cross_BH/F");
				}

			if (beampolartype==1) {
				SIM_Tree->Branch("cross_BH_x",&cross_BH_x ,"cross_BH_x/F");
				SIM_Tree->Branch("cross_BH_y",&cross_BH_y ,"cross_BH_y/F");
				}

			SIM_Tree->Branch("W_tot_unpol",&W_tot_unpol ,"W_tot_unpol/F");
			SIM_Tree->Branch("W_tot_pol", &W_tot_pol, "W_tot_pol/F");
			SIM_Tree->Branch("W_tot_pol_beam", &W_tot_pol_beam, "W_tot_pol_beam/F");
			SIM_Tree->Branch("W_tot_pol_target", &W_tot_pol_target, "W_tot_pol_target/F");
			SIM_Tree->Branch("W_TCS",&W_TCS ,"W_TCS/F");
			SIM_Tree->Branch("W_BH",&W_BH ,"W_BH/F");
			SIM_Tree->Branch("BSA",&BSA,"BSA/F");
			SIM_Tree->Branch("TSA",&TSA,"TSA/F");
			SIM_Tree->Branch("BTSA",&BTSA,"BTSA/F");
			SIM_Tree->Branch("target_spindir", &target_spindir,"target_spindir/I");
			SIM_Tree->Branch("beam_spindir", &beam_spindir,"beam_spindir/I");
			SIM_Tree->Branch("poltargetdeg", &poltargetdeg, "poltargetdeg/D");
			}

		SIM_Tree->Branch("polbeamdeg", &polbeamdeg, "polbeamdeg/D");
		SIM_Tree->Branch("thetamin_nocut", &thetamin_nocut, "thetamin_nocut/F");
		SIM_Tree->Branch("FlagSing",&FlagSing,"FlagSing/I");
		SIM_Tree->Branch("rad_event",&rad_event ,"rad_event/I");
		SIM_Tree->Branch("flux_rc",&flux_rc ,"flux_rc/F");
		SIM_Tree->Branch("EventNumber", &EventNumber , "EventNumber/L");
		SIM_Tree->Branch("TrueEventNumber", &TrueEventNumber , "TrueEventNumber/L");
		}

	/////////////////////////////////////////////////////////////////////////

	// Save user parameters
	for (int i=0;i<30;i++){
		param_initfile[i]=(float) param_init[i];
		}
	if (HEP==1 || HEP ==2){
		hepfile.open(Form("tcs_bh_output_hep_%d.dat", indexrun));//, std::ofstream::out);
		kinfile.open(Form("tcs_bh_output_kin_%d.dat", indexrun));}//, std::ofstream::out);
	else if (HEP==3 || HEP==4 || HEP==5 || HEP==6){
		hepfile.open(Form("tcs_bh_output_hep_%d.dat", indexrun));//, std::ofstream::out);
		}

	targetlenght_rc = targetlenght*0.5*rad_ob.UNX(Atarget, Ztarget);
	//targetlenght_rc = 0.0056; //targetlenght*0.5*rad_ob.UNX(Atarget, Ztarget);
	bremgam_integral = rad_ob.IntegralNph(Eelectron,Eelectron*0.0001, E_cutoff, targetlenght_rc);

	// tree with more info for debug
	//TreeInit(Dump_Tree);//, 0, 0);
	Dump_Tree->Branch("indexrun",&indexrun, "indexrun/I");
	Dump_Tree->Branch("param_initfile", &param_initfile,"param_initfile[30]/F");
	Dump_Tree->Branch("phase_space", &phase_space, "phase_space/F");
	Dump_Tree->Branch("TrueEventNumber", &TrueEventNumber , "TrueEventNumber/L");

	printf("\n Start the generation of events... \n");
	FirstEvent = 0;
	EventNumber = FirstEvent;

	printf("\n First event number: %lld",FirstEvent);
	TrueEventNumber=0;
	polarbeaminput = polbeamdeg;

	//********************* start generating event by event *********************************//
    for (int evv=0;evv<NTotEvents;evv++){
		TrueEventNumber++;

		// Initialisation of events
		VirtualFlag=-1; FlagSing=1; target_spindir=0; beam_spindir=0; check=0; i1=0; i2=0; i3=0;
		W_BH=0.; W_TCS=0.; W_tot_unpol=0.; W_tot_pol =0; W_tot_pol_beam = 0.; cross_tot_pol=0.; cross_tot_unpol=0; cross_tot_pol_beam=0; W_tot_pol_target=0; cross_tot_pol_target=0; cross_TCS=0;cross_BH=0; BSA=0.;  TSA=0.; BTSA=0.; cross_BH_x=0.; cross_BH_y=0; //W_DDVCS=0; W_VCS_born=0; W_VCS_nb=0; W_DVCS=0;
		thetamin_nocut=90; Flux_qr=0.; Flux_bmr=0.; theta_beam=0.; ttmin=0.; Psi_s = 0.; phi_s = 0; theta_s = 0; CosThetagg=0;

		CosThetaggCM=0; WW= 0; rad_event=0; d_internal=1;
		polpp=0.; polpm=0.; polmm=0.; polmp=0.; phisbin=-1; psisbin=-1; Ebin=-1; mtbin=-1;Qp2bin=-1; thetabin=-1; phibin=-1;//Q2bin=-1; Xbjbin=-1; PhiLHbin=-1;
		theta_gamma=0;   Gamma=0.; Beta=0.; PVirtual_CMeP=0; EVirtual_CMeP=0; test=0.;

		LV_minus_lab.SetPxPyPzE(0.,0.,0.,0.); LV_plus_lab.SetPxPyPzE(0.,0.,0.,0.); LV_gamma_lab.SetPxPyPzE(0.,0.,0.,0.); LV_Recoil_lab.SetPxPyPzE(0.,0.,0.,0.); LV_el_in.SetPxPyPzE(0.,0.,0.,0.); LV_el_out.SetPxPyPzE(0.,0.,0.,0.);  LV_plus.SetPxPyPzE(0.,0.,0.,0.); LV_minus.SetPxPyPzE(0.,0.,0.,0.);
		LV_Recoil_CMeP.SetPxPyPzE(0.,0.,0.,0.);  LV_Recoil.SetPxPyPzE(0.,0.,0.,0.); LV_target_lab.SetPxPyPzE(0.,0.,0.,0.938);// proton at rest

		for (int i=0;i<4;i++){
			ALV_minus_lab[i]=0.; ALV_plus_lab[i]=0.; ALV_gamma_in[i]=0.; ALV_el_in[i]=0.; ALV_Recoil_lab[i]=0.; ALV_el_out[i]=0.; //	ALV_gamma_out_lab[i]=0;
			}

		if (reaction==1 || reaction==11 ) check=RandomGen_TCS();

		if (check!=1){
			evv-=1;
			continue;
			}

		/////////////////////////////////////////////////////////////////////////
		// init beam options for TCS
		if (reaction==1 || reaction==11){
			polbeamdeg =pol.poltrans_elg(polarbeaminput, yy, beampolartype);
			if (beamtype==0){ // real photon beam
				Flux_qr=0; Flux_bmr=1; VirtualFlag=0; Q2=0;
				LV_gamma_lab.SetPxPyPzE(0,0,Egamma,Egamma);}
			else if (beamtype==1){ // electron beam and quasi-real photon
				Flux_qr= rad_ob.VFluxFactor(Q2max, Eelectron, Egamma);
				Flux_bmr= rad_ob.BremstrahlungSpectraNgammadiff(Eelectron,Egamma,Ztarget,Atarget,targetlenght);
				test =(rand() /(double)RAND_MAX);
				if (test < (Flux_qr/(Flux_qr+Flux_bmr)) ) {
					VirtualFlag=1;}
				else {
					VirtualFlag=0;
					Q2 = 0;
					}
				}
			}


		if (beamtype==1){
			LV_el_in.SetPxPyPzE(0.,0.,Eelectron,Eelectron);
			if (radcor>0){		// electron beam straggling + internal rad cor (real)
				LV_el_in = rad_ob.ElectronRC(LV_el_in,E_cutoff,Ztarget, targetlenght_rc,bremgam_integral, rad_event);
				if (radcor>1){
					d_internal = rad_ob.EqRad_lenght(Q2)*0.5;
					bremgam_integral_internal = rad_ob.IntegralNph(LV_el_in.E(),LV_el_in.E()*0.0001, E_cutoff, d_internal);
					LV_el_in = rad_ob.ElectronRC(LV_el_in,E_cutoff,1, d_internal,bremgam_integral_internal, rad_event);
					}
				}

			Ebeam = LV_el_in.E();
			if (reaction==32){ // recalculate Q2 in case of change in beam energy
				Q2 = 4.*Ebeam*(  Ebeam/ (1. + 2.*(Ebeam/M_Nucleon)*pow(sin(theta_beam*0.5),2)) )*pow(sin(theta_beam*0.5),2);
				}
			if (rad_event>0) {	// approximate correction of cross section vs E': ratio f_T(E')/f_T(E)
//				flux_rc = dis_ob.f_T_yQ(Egamma/Ebeam,Q2)/dis_ob.f_T_yQ(Egamma/Eelectron,Q2)*pow(Eelectron,2)/pow(Ebeam,2);
				// voir methode (ICI)
				}
			else flux_rc=1;


			if (VirtualFlag==1){
				theta_beam=2.*asin(sqrt( Q2/(4*Ebeam*(Ebeam-Egamma))));
				theta_gamma=2*PI-acos((Egamma+Q2/(2*Ebeam))/(pow((pow(Egamma,2)+Q2),0.5)));
				if (thetag_max>0.0001 && thetag_max*PI/180.<theta_gamma){
					evv -= 1;
					continue;
					}
				LV_gamma_lab.SetPxPyPzE(sqrt(pow(Egamma,2)+Q2)*sin(theta_gamma)*cos(phi_beam),
					sqrt(pow(Egamma,2)+Q2)*sin(theta_gamma)*sin(phi_beam), sqrt(pow(Egamma,2)+Q2)*cos(theta_gamma),Egamma);
				LV_el_out.SetPxPyPzE( (Ebeam-Egamma)*sin(theta_beam)*cos(phi_beam), (Ebeam-Egamma)*sin(theta_beam)*sin(phi_beam), (Ebeam-Egamma)*cos(theta_beam), Ebeam-Egamma);
				}
			else { // bremsstrahlung angle
				for (int tht=0;tht<8;tht++){
					theta_gamma = rad_ob.thetag_brem(Eelectron, Eelectron-Egamma,Ztarget);
					if (theta_gamma!=0) break;
					}
				if (thetag_max>0.0001 && thetag_max*PI/180.<theta_gamma){
					evv -= 1;
					continue;
					}

				theta_beam = PI-theta_gamma;
				LV_gamma_lab.SetPxPyPzE(sqrt(pow(Egamma,2)+Q2)*cos(phi_beam)*sin(theta_gamma),
					sqrt(pow(Egamma,2)+Q2)*sin(theta_gamma)*sin(phi_beam), sqrt(pow(Egamma,2)+Q2)*cos(theta_gamma), Egamma);
				LV_el_out.SetPxPyPzE(-LV_gamma_lab.Px(),-LV_gamma_lab.Py(),	Eelectron-Egamma*cos(theta_gamma), (float) Eelectron-Egamma);
				}


			epsilon = pol.fepsilon(Q2,Egamma,theta_beam);
			if (reaction==31 || reaction==32){
				LV_Recoil_lab.SetPxPyPzE( -LV_el_out.Px() ,-LV_el_out.Py(), Ebeam-LV_el_out.Pz(),Ebeam-LV_el_out.E()+M_Nucleon );// only non rad
				}

			// outgoing lepton straggling + internal real corrections
			if (radcor>0){
				bremgam_integral_internal = rad_ob.IntegralNph(LV_el_out.E(),LV_el_out.E()*0.0001, LV_el_out.E()*0.8, targetlenght_rc);
				LV_el_out = rad_ob.ElectronRC(LV_el_out,  LV_el_out.E()*0.8 ,Ztarget, targetlenght_rc,bremgam_integral_internal, rad_event);
				if (radcor>1){
					bremgam_integral_internal = rad_ob.IntegralNph(LV_el_out.E(),LV_el_out.E()*0.0001,  LV_el_out.E()*0.8, d_internal);
					LV_el_out = rad_ob.ElectronRC(LV_el_out,  LV_el_out.E()*0.8 ,1, d_internal,bremgam_integral_internal, rad_event);
					}
				}
			}


			//.......................... generate particle spin...............................................//
		if (reaction <=10 ) beam_spindir = pol.spin();
		if (reaction==1 ){ // tcs
			target_spindir=pol.spin();
			phi_s = pol.fphis(phi_beam,targetpoldir);
			theta_s = pol.fthetas(theta_gamma, targetpoldir);
			if ((reaction==1 || reaction==8) && beampolartype>0){
				Psi_s = pol.fpsis(phi_beam, beampolartype);
				}
			}

		// generate 4-vectors lab->CM->decay->CM->lab->axis rotations (real lab)

		if (reaction < 30){
			test = CMlabT.get_gammaPout(LV_gamma_out_lab, LV_Recoil_lab, Beta, Gamma, WW, CosThetaggCM, EVirtual_CMeP, PVirtual_CMeP,
				(float) Ebeam, (float) Egamma,(float) Q2, (float) Qp2, (float) tt, (float) thetagg_min, (float) thetagg_max);

			if (test==0) { evv-=1; continue; } // warning! if bug: continuous loop

			// Final lepton pair production from virtual photon rest frame
			if (reaction==1 || reaction==11 ){
				pair.get_pair(LV_minus_lab,LV_plus_lab, EVirtual_CMeP, PVirtual_CMeP, acos(CosThetaggCM), Theta_CMV, Phi_CMV, M_lepton);	
				part_op.BoostBack(LV_minus_lab, Beta, Gamma);
				part_op.BoostBack(LV_plus_lab, Beta, Gamma);
				}

			// electroproduction: physics angle, go back to 'real' lab. All: random rotation reaction plane
			if ((reaction==1 && beamtype ==1) || (reaction==11 && beamtype==1)){
				CMlabT.rotate_outelectro_backlab(LV_Recoil_lab, phi_beam, theta_gamma, Phi_LH);
				CMlabT.rotate_outelectro_backlab(LV_plus_lab, phi_beam, theta_gamma, Phi_LH);
				CMlabT.rotate_outelectro_backlab(LV_minus_lab, phi_beam, theta_gamma, Phi_LH);
				}
			else  if ((reaction==1 && beamtype==0) || (reaction==11 && beamtype==0))  {
				part_op.Rot_clock_Z(LV_Recoil_lab,phi_beam);
				part_op.Rot_clock_Z(LV_minus_lab,phi_beam);
				part_op.Rot_clock_Z(LV_plus_lab,phi_beam);
				}
			if (radcor>0 && outlepton==1){
				bremgam_integral_internal = rad_ob.IntegralNph(LV_plus_lab.E(),LV_plus_lab.E()*0.0001,LV_plus_lab.E()*0.8, targetlenght_rc);
				LV_plus_lab = rad_ob.ElectronRC(LV_plus_lab ,LV_plus_lab.E()*0.8 ,Ztarget, targetlenght_rc,bremgam_integral_internal, rad_event);
				bremgam_integral_internal = rad_ob.IntegralNph(LV_minus_lab.E(),LV_minus_lab.E()*0.0001,LV_minus_lab.E()*0.8, targetlenght_rc);
				LV_minus_lab = rad_ob.ElectronRC(LV_minus_lab ,LV_minus_lab.E()*0.8 ,Ztarget, targetlenght_rc,bremgam_integral_internal, rad_event);
				if (radcor>1){
					d_internal = rad_ob.EqRad_lenght(Qp2)*0.5;
					bremgam_integral_internal = rad_ob.IntegralNph(LV_plus_lab.E(),LV_plus_lab.E()*0.0001,LV_plus_lab.E()*0.8, targetlenght_rc);
					LV_plus_lab = rad_ob.ElectronRC(LV_plus_lab ,LV_plus_lab.E()*0.8 ,1, d_internal,bremgam_integral_internal, rad_event);
					bremgam_integral_internal = rad_ob.IntegralNph(LV_minus_lab.E(),LV_minus_lab.E()*0.0001,LV_minus_lab.E()*0.8, targetlenght_rc);
					LV_minus_lab = rad_ob.ElectronRC(LV_minus_lab ,LV_minus_lab.E()*0.8 ,1, d_internal,bremgam_integral_internal, rad_event);
					}
				}

				// fill arrays for output
			part_op.FillArray_LV(LV_Recoil_lab,ALV_Recoil_lab);
			part_op.FillArray_LV(LV_gamma_lab,ALV_gamma_in);
			if (beamtype!=0){
				part_op.FillArray_LV(LV_el_out, ALV_el_out);
				part_op.FillArray_LV(LV_el_in, ALV_el_in);
				}
			if (reaction==1 || reaction==11 ){
				part_op.FillArray_LV(LV_minus_lab, ALV_minus_lab);
				part_op.FillArray_LV(LV_plus_lab, ALV_plus_lab);
				}
			}	// end of filling 4 vectors and arrays


		// lab frame cuts to apply here (speed/size optimization)
		if (HEP==5 || HEP==6){
			if (reaction==1 || reaction==11 ){
				if (mom_shms>0.01 && (LV_minus_lab.P()>mom_shms*1.25 || LV_minus_lab.P()<mom_shms*0.75)) { evv-=1; continue; }
				if (mom_hms>0.01 && (LV_plus_lab.P()>mom_hms*1.25 || LV_plus_lab.P()<mom_hms*0.75)){ evv -=1; continue; }
				if (theta_shms>0.01 && (LV_minus_lab.Theta()>(theta_shms+3)*PI/180. || LV_minus_lab.Theta()<(theta_shms-3)*PI/180.)){ evv-=1; continue; }
				if (theta_hms>0.01 && (LV_plus_lab.Theta()>(theta_hms+3)*PI/180. || LV_plus_lab.Theta()<(theta_hms-3)*PI/180.)){ evv-=1; continue; }
				}
			}

		EventNumber++; // pure phase-space end
		if (reaction==1 || reaction==11){ 	// flag for phase space cut in TCS
			i1=(int)SetBins( Egamma,thmEb,thmE);
			i2=(int)SetBins( (-tt) ,thmTb,thmT);
			i3=(int)SetBins( Qp2,thmQb,thmQ);
				if (Egamma<thmEb[0]) {i1=0;} if (Egamma>thmEb[thmE-1]) {i1=thmE-1;}
				if (-tt<thmTb[0]) {i2=0;} if (-tt>thmTb[thmT-1]) {i2=thmT-1;}
				if (Qp2<thmQb[0]) {i3=0;} if (Qp2>thmQb[thmQ-1]) {i3=thmQ-1;}
				if (!(i1<0 || i1>=thmE || i2<0 || i2>=thmT || i3<0 || i3>=thmQ)){
					thetamin_nocut = thm[(int) i1][(int) i2][(int) i3];
					}
				if ((Theta_CMV*180./(float) PI)<thetamin_nocut || (Theta_CMV*180./(float) PI)>(180.-thetamin_nocut) ){
					if (Theta_CMV<PI/2.){
						if (Phi_CMV<2.62 || Phi_CMV>3.67)  FlagSing=0;
						else FlagSing=1;
						}
					else {
						if (Phi_CMV<0.52 || Phi_CMV>5.76) FlagSing=1;
						else FlagSing=0;
						}
					}
				else FlagSing=0;
		}

		// event weighting
		if (reaction==1){
			Ebin = (int) SetBins(Egamma,Egamma_bintable,NEB);
			mtbin = (int) SetBins(-tt,tt_bintable,NT);
			Qp2bin = (int) SetBins(Qp2,Qp2_bintable,NQp2);
			thetabin = (int) SetBins(Theta_CMV,theta_bintable,NTh);
			phibin = (int) SetBins(Phi_CMV,phi_bintable,NPhi);
			if (targetpoldir==1 || targetpoldir==2) phisbin = (int) SetBins(phi_s,phis_bintable,NPhis);
			else phisbin=0;
			if (beampolartype==1) psisbin = (int)SetBins(Psi_s,psis_bintable,NPsis);
			else psisbin=0;

			// check and modify binning (ICI)
			// read cross section tables

			if (targetpoldir==1 || targetpoldir==2){
				if (phi_s<PI) {
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polpp, polpm, polmp, polmm, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin,(int) phisbin );
					}
				else {
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polpm, polpp, polmm, polmp, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin,(int) (phisbin-NPhis/2) );
					}
				}
			else if (beampolartype==1){
				if (Psi_s<PI){
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polpp, polpm, polmp, polmm, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin, (int) psisbin );
					}
				else if (Psi_s>=PI/2. && Psi_s<PI){
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polmp, polmm, polpp, polpm, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin, (int) (psisbin-NPsis/4) );
					}
				else if (Psi_s>=PI && Psi_s<0.75*PI){
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polpp, polpm, polmp, polmm, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin, (int) (psisbin-NPsis/2) );
					}
				else if (Psi_s>=PI){
					cross_tot_unpol= linear_interpol_tcs6D_4( ddvcs_table_provider, cross_BH, cross_TCS, polmp, polmm, polpp, polpm, cross_BH_x, cross_BH_y, -tt, Qp2, Phi_CMV, Theta_CMV, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, targetpoldir, beampolartype, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin, (int) (psisbin-(3*NPsis/4)) );
					}
				}
			else {
				cross_tot_unpol= linear_interpol_tcs5D( ddvcs_table_provider, cross_BH, cross_TCS, polpp, polpm, polmp, polmm, Egamma, -tt, Qp2, Phi_CMV, Theta_CMV, Egamma_bintable, tt_bintable, Qp2_bintable, theta_bintable, phi_bintable, (int) Ebin, (int) mtbin, (int) Qp2bin, (int) thetabin, (int) phibin );
				}
			if (polpp<=0 || polpm<=0 || polmp<=0 || polmm<=0) {
				evv-=1; // case of out the range after interpolation (very rare), binning effect
				continue;
				}

			// polarized cross section with dilution factor and spin orientation
			cross_tot_pol_beam = pol.cross_poldilut( (polpp + polpm)*0.5, (polmp+polmm)*0.5,beam_spindir,polbeamdeg);
			if (targetpoldir>0) {
				cross_tot_pol_target = pol.cross_poldilut( (polpp + polmp)*0.5, (polmm+polpm)*0.5,target_spindir,poltargetdeg);
				cross_tot_pol = pol.cross_doublepoldilut(polpp, polpm, polmp, polmm, beam_spindir, target_spindir, polbeamdeg, poltargetdeg);
				}
			else {
				cross_tot_pol_target = cross_tot_unpol;
				cross_tot_pol=cross_tot_pol_beam;
				}

			// conversion to picobarn and phase space in dtheta
			cross_tot_unpol *=1000*sin(Theta_CMV);
			cross_tot_pol *=1000*sin(Theta_CMV);
			cross_tot_pol_beam *=1000*sin(Theta_CMV);
			cross_tot_pol_target *=1000*sin(Theta_CMV);
			cross_TCS *=1000*sin(Theta_CMV);
			cross_BH *= 1000*sin(Theta_CMV);
			cross_BH_x*=1000*sin(Theta_CMV); cross_BH_y*=1000*sin(Theta_CMV);

			// approximate asymmetries for specific event kinematic, assume no dilution
			BSA = (polpp+polpm-polmp-polmm)/(polpp+polpm+polmp+polmm);
			TSA = (polpp+polmp-polpm-polmm)/(polpp+polpm+polmp+polmm);
			BTSA = (polpp+polmm-polmp-polpm)/(polpp+polpm+polmp+polmm);

			// weights accounting e-gamma coupling
			W_tot_pol = cross_tot_pol * (Flux_qr+Flux_bmr);
			W_tot_pol_beam = cross_tot_pol_beam * (Flux_qr+Flux_bmr);
			W_tot_pol_target = cross_tot_pol_target * (Flux_qr+Flux_bmr);
			W_tot_unpol=cross_tot_unpol*(Flux_qr+Flux_bmr);
			W_BH=cross_BH*(Flux_qr+Flux_bmr);
			W_TCS=cross_TCS*(Flux_qr+Flux_bmr);

			if (HEP==1 || HEP==2){
				if (beamtype==0){
					hepfile << "3" << endl;
					hepfile << "1 " << "11" << " 0 0 " << ALV_minus_lab[1] << " "<< ALV_minus_lab[2] << " " << ALV_minus_lab[3] << " "<< m_el << endl;
					hepfile << "1 " << "-11" << " 0 0 "<< ALV_plus_lab[1] << " " << ALV_plus_lab[2] << " " << ALV_plus_lab[3] << " "<< m_el << endl;
					hepfile << "1 " << "2212" << " 0 0 " << ALV_Recoil_lab[1] << " "<< ALV_Recoil_lab[2] << " " << ALV_Recoil_lab[3] << " " << M_Nucleon << endl;
					kinfile<< Qp2 << " " << tt << " " << ttmin << " " << Egamma <<" " << Phi_CMV << " " << Theta_CMV << " " <<CosThetagg<<" "<<phi_beam<< " " <<W_tot_unpol << " " <<W_tot_pol<<" "<<W_tot_pol_beam<<" "<<W_tot_pol_target<<" "<<W_TCS<<" "<<W_BH<<" "<<BSA<<" "<<TSA<<" "<<BTSA<<" "<<target_spindir<<" "<<beam_spindir<<" "<<poltargetdeg<<" "<<polbeamdeg<<" " <<thetamin_nocut<<" "<<FlagSing<<" "<<rad_event<<" "<<TrueEventNumber;
					if (beampolartype==0){
						if (targetpoldir==0 || targetpoldir==3) kinfile<<" "<<endl;
						else kinfile<<" "<<phi_s<<endl;
						}
					else {
						kinfile<<" "<<cross_BH_x<<" "<<cross_BH_y<<Psi_s<<endl;
						}
					}
				else {
					hepfile << "4" << endl;
					hepfile << "1 " << "11" << " 0 0 " << ALV_minus_lab[1] << " "<< ALV_minus_lab[2] << " " << ALV_minus_lab[3] << " "<< m_el << endl;
					hepfile << "1 " << "-11" << " 0 0 "<< ALV_plus_lab[1] << " " << ALV_plus_lab[2] << " " << ALV_plus_lab[3] << " "<< m_el << endl;
					hepfile << "1 " << "2212" << " 0 0 " << ALV_Recoil_lab[1] << " "<< ALV_Recoil_lab[2] << " " << ALV_Recoil_lab[3] << " " << M_Nucleon << endl;
					hepfile << "1 " << "11" << " 0 0 " << ALV_el_out[1] << " "<< ALV_el_out[2] << " " << ALV_el_out[3] << " " << m_el << endl;
					kinfile<< Qp2 << " " << tt << " " << ttmin << " " << Egamma <<" " << Phi_CMV << " " << Theta_CMV << " " <<CosThetagg<<" "<<phi_beam<< " "
					<< Q2<<" "<<yy<<" "<<epsilon<<" "<<VirtualFlag<<" "<<Flux_qr<<" "<<Flux_bmr<<" "<<theta_gamma<<" "<<theta_beam
					<<" "<<cross_tot_unpol << " " <<cross_tot_pol<<" "<<cross_tot_pol_beam<<" "<<cross_tot_pol_target<<" "<<cross_TCS<<" "<<cross_BH
					<<" "<<W_tot_unpol << " " <<W_tot_pol<<" "<<W_tot_pol_beam<<" "<<W_tot_pol_target<<" "<<W_TCS<<" "<<W_BH<<" "<<BSA<<" "<<TSA<<" "<<BTSA<<" "<<target_spindir<<" "<<beam_spindir<<" "<<poltargetdeg<<" "<<polbeamdeg<<" " <<thetamin_nocut<<" "<<FlagSing<<" "<<rad_event<<" "<<TrueEventNumber;
					if (beampolartype==0){
						if (targetpoldir==0 || targetpoldir==3) kinfile<<" "<<endl;
						else kinfile<<" "<<phi_s<<endl;
						}
					else {
						kinfile<<" "<<cross_BH_x<<" "<<cross_BH_y<<Psi_s<<endl;
						}
					}
				}
			else if (HEP==3 || HEP==4 || HEP==5 || HEP==6){
				hepfile << -ALV_plus_lab[2] << " "<< ALV_plus_lab[1] << " " << ALV_plus_lab[3] << " " << ALV_plus_lab[0] << " 0 "<< -ALV_minus_lab[2] << " " << ALV_minus_lab[1] << " " << ALV_minus_lab[3] << " "<< ALV_minus_lab[0] << " 0 " << W_BH << " "<<
				W_TCS <<" "<<W_tot_unpol<<" "<<BSA<<" "<<TSA<<" "<<BTSA<<" "<<(Flux_qr+Flux_bmr)<<" "<<Ebeam<<endl;
				}
			}

		else if (reaction==11 ){
			// no weighting, just phase space for TCS
			if (HEP==1 || HEP==2){
				hepfile << "3" << endl;
				hepfile << "1 " << "11" << " 0 0 " << ALV_minus_lab[1] << " "<< ALV_minus_lab[2] << " " << ALV_minus_lab[3] << " "<< m_el << endl;
				hepfile << "1 " << "-11" << " 0 0 "<< ALV_plus_lab[1] << " " << ALV_plus_lab[2] << " " << ALV_plus_lab[3] << " "<< m_el << endl;
				hepfile << "1 " << "2212" << " 0 0 " << ALV_Recoil_lab[1] << " "<< ALV_Recoil_lab[2] << " " << ALV_Recoil_lab[3] << " " << M_Nucleon << endl;
				kinfile<< Qp2 << " " << tt << " " << ttmin << " " << Egamma <<" " << Phi_CMV << " " << Theta_CMV << " " <<CosThetagg<<" "<<phi_beam<< " "<<thetamin_nocut<<" "<<FlagSing<<" "<<rad_event<<endl;
				}
			else if (HEP==3 || HEP==4 || HEP==5 || HEP==6){
				hepfile << -ALV_plus_lab[2] << " "<< ALV_plus_lab[1] << " " << ALV_plus_lab[3] << " " << ALV_plus_lab[0] << " 0 "<< -ALV_minus_lab[2] << " " << ALV_minus_lab[1] << " " << ALV_minus_lab[3] << " "<< ALV_minus_lab[0] << " 0 " <<"1" << " "<<"1" <<" "<<"1"<<" "<<-ALV_Recoil_lab[2]<<" "<<ALV_Recoil_lab[1]<<" "<<ALV_Recoil_lab[3]<<" "<<flux_rc<<" "<<Ebeam<<endl;
				}
			}

		// put back angles from 0 to 2pi
		if (Phi_LH>2.*PI) Phi_LH -= 2.*PI;
		if (Phi_CMV>2.*PI) Phi_CMV -= 2.*PI;
		//if (phis>2.*PI) phis -= 2.*PI;

		if (HEP==0 || HEP==1 || HEP==4 || HEP==6) SIM_Tree->Fill();// fill the root file
		if (evv==NTotEvents-1) Dump_Tree->Fill();	// keep 50 first events for double checks and studies

		//HDDM STUFF
		tmpEvt_t tmpEvt;
		tmpEvt.beam = LV_gamma_lab;
		tmpEvt.target = LV_target_lab;
		tmpEvt.q1 = LV_minus_lab;
		tmpEvt.q2 = LV_plus_lab;
		tmpEvt.recoil = LV_Recoil_lab;
		tmpEvt.nGen = 3;
		tmpEvt.rxn = 2;  //the rxn number 2 is for proton as recoil in HDDM i.e HddmOut.h
		//tmpEvt.weight = fullWeight;
		hddmGo.write(tmpEvt,runNum,TrueEventNumber);

		if (evv%1000==0) cout<<" "<<evv<<" events ... ";

		}
		//********************* end of generating event by event *********************************//

	logfile<<indexrun<<" "<<NTotEvents<<" "<<TrueEventNumber<<" "<<phase_space<<endl;
	logfile<<"Above: run_index, N_tot_file, N_ActualGen, PhaseSpace:"<<endl;
	logfile<<"Below: gen param (order of input file):"<<endl;

	for (int i=0;i<30;i++){
		logfile<<param_init[i]<<" ";
		}
	logfile<<"process: "<<reaction<<" Input files version "<<endl;
    string ff; ff = Form("set.csh");
	ifstream inset; inset.open(ff.c_str());
	if (inset){
		while (inset.good()){
			getline(inset , ff);
			logfile<<ff<<endl;
			}
		}

	cout<<"Write output file"<<endl;
	file->Write();
	file->Close();

	cout<<"Generation is done"<<endl;
	cout<<"Output file is tcs_bh_output_"<<indexrun<<".root"<<endl;
	cout<<"SIM_Tree = generated events"<<endl;
	cout<<"Dump_Tree = last event info"<<endl;

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Compute duration
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

	logfile<<"Elapsed time: "<< elapsed.count() << " seconds" <<endl;

	return 100;
} // end main
