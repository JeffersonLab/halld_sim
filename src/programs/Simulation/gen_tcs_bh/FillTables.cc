#define FillTables_cxx
#include "FillTables.h"
#include <limits>
#include <unistd.h>
using namespace std;

int FillTables::FillTable(int reaction, TableProvider& table_provider){

	genSettings_t genSettings0;
	TString Tablepath0 = genSettings0.xsecTablepath;
        int Qp2_LimitType0 = genSettings0.Qp2_LimitType;
        printf("Start to fill all the tables ...");

        ifstream infile1x,infile1z;
        string ff1x, ff1z, spoub;
        char proc[24];
        double em=0.,tm=0.,Qpm=0.,Thm=0.,phim=0.,psim=0.; //Qm=0.,PLm=0.,
        int ii=0, jj=0, kk=0, ll=0, mm=0,nn=0,countercheck=0; //oo=0,
        int nn_xx=-1,nn_tt=-1,nn_qp=-1,nn_th=-1,nn_pc=-1,nn_ps=-1; //,nn_qq=-1,nn_pl=-1
        double res[15]={0.};
        float CosThetaMax=1.;

        if (reaction==1) strcpy(proc, "tcs");

        cout<<"reaction index: "<<reaction<<endl;

	if (strcmp(proc, "tcs")==0){

                if (protonorneutron==1){
                        if (targetpoldir==1 || targetpoldir==2){
                                cout<<"Read table for circularly polarized photon and transversely polarized target"<<endl;
                                //ff1x=Form(Tablepath0+"/grid_table_tcs_circperp.dat");
                                if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_circperp_LQ2.dat");}
                                if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_circperp_HQ2.dat");}
                                if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_circperp_FQ2.dat");}
                                //ff1x=Form("Data/grid_table_tcs_circperp.dat");
                        } else {
                                if (beampolartype==1){
                                        cout<<"Read table for linearly polarized photon and longitudinally polarized target"<<endl;
                                        //ff1x=Form(Tablepath0+"/grid_table_tcs_linlong.dat");
                                        if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_linlong_LQ2.dat");}
                                        if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_linlong_HQ2.dat");}
                                        if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_linlong_FQ2.dat");}
                                        //ff1x=Form("Data/grid_table_tcs_linlong.dat");
                                } else {
                                        cout<<"Read table for circularly polarized photon and longitudinally polarized target"<<endl;
                                        //ff1x=Form(Tablepath0+"/grid_table_tcs_circlong.dat");
                                        if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_circlong_LQ2.dat");}
                                        if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_circlong_HQ2.dat");}
                                        if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_circlong_FQ2.dat");}
                                        //ff1x=Form("Data/grid_table_tcs_circlong.dat");
                                }
                        }
                } else if (protonorneutron==2){
                        cout<<"Read cross section off neutron"<<endl;
                        if (targetpoldir==1 || targetpoldir==2){
                                cout<<"Read table for circularly polarized photon and transversely polarized target"<<endl;
                                //ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circperp.dat");
                                if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circperp_LQ2.dat");}
                                if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circperp_HQ2.dat");}
                                if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circperp_FQ2.dat");}
                                //ff1x=Form("Data/grid_table_tcs_neutron_circperp.dat");
                        } else {
                                if (beampolartype==1){
                                        cout<<"Read table for linearly polarized photon and longitudinally polarized target"<<endl;
                                        //ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_linlong.dat");
                                        if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_linlong_LQ2.dat");}
                                        if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_linlong_HQ2.dat");}
                                        if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_linlong_FQ2.dat");}
                                        //ff1x=Form("Data/grid_table_tcs_neutron_linlong.dat");
                                } else {
                                        cout<<"Read table for circularly polarized photon and longitudinally polarized target"<<endl;
                                        //ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circlong.dat");
                                        if (Qp2_LimitType0==1){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circlong_LQ2.dat");}
                                        if (Qp2_LimitType0==2){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circlong_HQ2.dat");}
                                        if (Qp2_LimitType0==3){ff1x=Form(Tablepath0+"/grid_table_tcs_neutron_circlong_FQ2.dat");}
                                        //ff1x=Form("Data/grid_table_tcs_neutron_circlong.dat");
                                }
                        }

                } else {
                        cout<<"ERROR: target = proton or neutron"<<endl;
                        return 5;
                }
                infile1x.open(ff1x.c_str());
                if (!infile1x) {
                        cout<<"no input for TCS cross section"<<endl;
                        cout<<"EXIT program with ERROR"<<endl;
                        return 0;
                }
                if (targetpoldir==1 || targetpoldir==2){
                        infile1x>>nn_xx>>nn_tt>>nn_qp>>nn_th>>nn_pc>>nn_ps>>CosThetaMax;
                } else {
                        if (beampolartype==1){
                                infile1x>>nn_xx>>nn_tt>>nn_qp>>nn_th>>nn_pc>>nn_ps>>CosThetaMax;
                        } else {
                                infile1x>>nn_xx>>nn_tt>>nn_qp>>nn_th>>nn_pc>>CosThetaMax;
                        }
                }

                if  (nn_xx!=NEB || nn_tt!=NT || nn_qp!=NQp2
                        || nn_th!=NTh || nn_pc!=NPhi){
                        cout<<"WARNING: not the same number of bins in table and in generator"<<endl;
                        cout<<"->check TreatOptions.h"<<endl;
                        cout<<"Eb: "<<nn_xx<<" "<<NEB<<endl;
                        cout<<"-t: "<<nn_tt<<" "<<NT<<endl;
                        cout<<"Qp: "<<nn_qp<<" "<<NQp2<<endl;
                        cout<<"Th: "<<nn_th<<" "<<NTh<<endl;
                        cout<<"Ph: "<<nn_pc<<" "<<NPhi<<endl;
                } else cout<<"initialization, N bins OK"<<endl;
                if (targetpoldir==1 || targetpoldir==2){
                        if (nn_ps != NPhis){
                                cout<<"WARNING: not the same number of bins in table and in generator"<<endl;
                                cout<<"phi_s: "<<nn_ps<<" "<<NPhis<<endl;
                        }
                }
                if (beampolartype==1){
                        if (nn_ps != NPsis){
                                cout<<"WARNING: not the same number of bins in table and in generator"<<endl;
                                cout<<"psi_s: "<<nn_ps<<" "<<NPsis<<endl;
                        }
                }
                if (beampolartype==1 && (targetpoldir==1 || targetpoldir==2)){
                        cout<<"WARNING: transversaly polarized target + linearly polarized beam is not included in TCS event generator"<<endl;
                        cout<<"photon beam is switch to circular polarization"<<endl;
                        beampolartype=0;
                }
                std::cout<<"Max theta_gamma-gamma = "<<CosThetaMax<<std::endl;
                float Egamma_bintable[NEB+1], tt_bintable[NT+1], Qp2_bintable[NQp2+1],phi_bintable[NPhi+1],theta_bintable[NTh+1], Phis_bintable[NPhis+1], Psis_bintable[NPsis+1];
                for (int i=0; i<=NEB; i++){
                      Egamma_bintable[i] = (float) LinearBins(NEB,EMINI,EMAXI,i);
                }
                for (int i=0; i<=NT; i++){
                      tt_bintable[i] = (float) LinearBins(NT,TMINI,TMAXI,i);
                }
                for (int i=0; i<=NQp2; i++){
                        Qp2_bintable[i] = (float) LinearBins(NQp2,QP2MINI,QP2MAXI,i);
                }
                for (int i=0; i<=NPhi; i++){
                        phi_bintable[i] = (float) LinearBins(NPhi,PHIMINI,PHIMAXI,i);
                }
                for (int i=0; i<=NTh; i++){
                        theta_bintable[i] = (float) LinearBins(NTh,THMINI,THMAXI,i);
                }
                if (targetpoldir==1 || targetpoldir==2){
                     for (int i=0; i<=NPhis; i++){
                        Phis_bintable[i] = (float) LinearBins(NPhis,PHISMINI,PHISMAXI,i);
                     }
                }
                if (beampolartype==1){
                     for (int i=0; i<=NPsis; i++){
                        Psis_bintable[i] = (float) LinearBins(NPsis,PSISMINI,PSISMAXI,i);
                     }
                }

                cout<<"\n start: insert values associating cross sections to kinematics"<<endl;
                cout<<">>> this step can take few minutes"<<endl;
                ii=0;jj=0;kk=0;ll=0;mm=0;nn=0;
                while (infile1x.good()){
                        if (!(infile1x>>em)) break;
                        if (targetpoldir==1 || targetpoldir==2){
                                infile1x>> tm>>Qpm>>Thm>>phim>>psim>> res[0]>>res[1]>>res[2]>>res[3]>>res[4]>>res[5];
                        } else {
                                if (beampolartype==1){
                                        infile1x>> tm>>Qpm>>Thm>>phim>>psim>> res[1]>>res[2]>>res[3]>>res[4]>>res[5]>>res[6]>>res[7];

                                } else {
                                        infile1x>> tm>>Qpm>>Thm>>phim>> res[0]>>res[1]>>res[2]>>res[3]>>res[4]>>res[5];
                                }
                        }

                        ii=SetBins(em+0.001,Egamma_bintable,NEB);
                        jj=SetBins(tm+0.001,tt_bintable,NT);
                        kk=SetBins(Qpm+0.001,Qp2_bintable,NQp2);
                        ll=SetBins(Thm+0.001,theta_bintable,NTh);
                        mm=SetBins(phim+0.001,phi_bintable,NPhi);
                        if (beampolartype==1) nn=SetBins(psim+0.001,Psis_bintable,NPsis);
                        if (targetpoldir==1 || targetpoldir==2) nn=SetBins(psim+0.001,Phis_bintable,NPhis);

                        if (targetpoldir==1 || targetpoldir==2 || beampolartype==1){

                            if(!table_provider.insert_tcs_pol( {ii,jj,kk,ll,mm,nn}, {res[0],res[1],res[2],res[3],res[4],res[5], res[6], res[7] })) {
                                                                std::cout << "[Fill Table TCS pol] FAILED TO INSERT " << Vars_Cin_tcs_pol({ii,jj,kk,ll,mm,nn})
                                                                << " -> " << Vars_Sef_tcs_pol({res[0],res[1],res[2],res[3],res[4],res[5]}) << std::endl;
return 5;
                            } else {
                                countercheck++;
                            }
                        } else {
                            if(!table_provider.insert_tcs( {ii,jj,kk,ll,mm}, {res[0],res[1],res[2],res[3],res[4],res[5]})) {
                                                                std::cout << "[Fill Table] FAILED TO INSERT " << Vars_Cin_tcs({ii,jj,kk,ll,mm})
                                                                << " -> " << Vars_Sef_tcs({res[0],res[1],res[2],res[3],res[4],res[5]}) << std::endl;
return 5;
                            } else {
                                countercheck++;
                            }
                        }

                }
                cout<<"size of table: "<<countercheck<<endl;
                infile1x.close();
                return 1;

	} // end tcs
	 else {
                printf ("\n process= use lower case call help");
                cout<<"EXIT with ERROR"<<endl;
                return 0;
        }

} // end fill table










