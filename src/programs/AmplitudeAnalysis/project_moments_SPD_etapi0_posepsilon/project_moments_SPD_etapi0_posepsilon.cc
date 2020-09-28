#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>

#include "IUAmpTools/FitResults.h"

#include "TFile.h"

using namespace std;






int main( int argc, char* argv[] ){
    
    // these params should probably come in on the command line
    double lowMass = 0.7;
    double highMass = 3.0;
    enum{ kNumBins = 45 };
    double lowt = 0;
    double hight = 1.2;
    enum{ kNumBinst = 4 };
    string fitDir("etaprimepi0");
                     
    // set default parameters
    
    string outfileName("");
    
    // parse command line
    
    for (int i = 1; i < argc; i++){
        
        string arg(argv[i]);
        
        if (arg == "-o"){
            if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
            else  outfileName = argv[++i]; }
        if (arg == "-h"){
            cout << endl << " Usage for: " << argv[0] << endl << endl;
            cout << "\t -o <file>\t Ouput text file" << endl;
            exit(1);}
    }
    
    if (outfileName.size() == 0){
        cout << "No output file specified" << endl;
        exit(1);
    }
    
    double step = ( highMass - lowMass ) / kNumBins;
    double stept = ( hight - lowt ) / kNumBinst;
    
    ofstream outfile;
    outfile.open( outfileName.c_str() );
    
    outfile <<"M"<<"\t"<<"t"
	<<"\t"<<"H1_00"<<"\t"<<"H1_00uncert."
        <<"\t"<<"H0_00"<<"\t"<<"H0_00uncert."
	<<"\t"<<"H1_10"<<"\t"<<"H1_10uncert."
        <<"\t"<<"H0_10"<<"\t"<<"H0_10uncert." 
	<<"\t"<<"H1_11"<<"\t"<<"H1_11uncert."
	<<"\t"<<"H0_11"<<"\t"<<"H0_11uncert."
	<<"\t"<<"H1_20"<<"\t"<<"H1_20uncert."
	<<"\t"<<"H0_20"<<"\t"<<"H0_20uncert."
	<<"\t"<<"H1_21"<<"\t"<<"H1_21uncert."
	<<"\t"<<"H0_21"<<"\t"<<"H0_21uncert."
	<<"\t"<<"H0_22"<<"\t"<<"H0_22uncert."
	<<"\t"<<"H1_22"<<"\t"<<"H1_22uncert."<<endl;



    // descend into the directory that contains the bins
    chdir( fitDir.c_str() );
    
    for( int i = 0; i < kNumBins;i++ ){
      for( int j = 0; j < kNumBinst; j++ ){  
	cout<<"bin "<<i<<"_"<<j<<endl;
        
        ostringstream dir;
        dir << "bin_" << i<<"_"<<j;
        chdir( dir.str().c_str() );
        
        ostringstream resultsFile;
        resultsFile << "bin_" << i <<"_"<<j<<".fit";
        
     

        FitResults results( resultsFile.str() );
        
        if( !results.valid() ){
	  outfile<< lowMass + step * i + step / 2. << "\t";
	  outfile<< lowt + stept * j + stept / 2. << "\t";
	  outfile <<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"
		  <<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"
		  <<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;            
            chdir( ".." );
            continue;
        }
        
        // print out the bin center
        outfile << lowMass + step * i + step / 2. << "\t";
	outfile << lowt + stept * j + stept / 2. << "\t";



	//complex<double>S0pRe_comp(0,0);
	complex<double>S0pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::S0+");
	complex<double>P0pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::P0+");
	complex<double>P1pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::P1+");
	complex<double>D0pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::D0+");
	complex<double>D1pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::D1+");
	complex<double>D2pRe_comp=results.scaledProductionParameter("EtaPrimePi0::PositiveRe::D2+");
                                    


	
	pair <double,double> H1_00(2.*(pow(abs(S0pRe_comp),2.)+pow(abs(P0pRe_comp),2.)+pow(abs(D0pRe_comp),2.)),0.);
	pair <double,double> H0_00(H1_00.first+2.*(pow(abs(P1pRe_comp),2.)+pow(abs(D1pRe_comp),2.)+pow(abs(D2pRe_comp),2.)),0.);
	pair <double,double> H1_10((8./sqrt(15.))*real(P0pRe_comp*conj(D0pRe_comp))+(4./sqrt(3.))*real(S0pRe_comp*conj(P0pRe_comp)),0.);
       	pair <double,double> H0_10(H1_10.first+(4./sqrt(5.))*real(P1pRe_comp*conj(D1pRe_comp)),0.);

        pair <double,double> H1_11((2./sqrt(5.))*real(P0pRe_comp*conj(D1pRe_comp))-(2./sqrt(15.))*real(P1pRe_comp*conj(D0pRe_comp))+(2./sqrt(3.))*real(S0pRe_comp*conj(P1pRe_comp)),0.);
        pair <double,double> H0_11(H1_11.first+2.*sqrt(2./5.)*real(P1pRe_comp*conj(D2pRe_comp)),0.);
        pair <double,double> H1_20((4./5.)*pow(abs(P0pRe_comp),2.)+(4./7.)*pow(abs(D0pRe_comp),2.)+(4./sqrt(5.))*real(S0pRe_comp*conj(D0pRe_comp)),0.);
	pair <double,double> H0_20(H1_20.first-(2./5.)*pow(abs(P1pRe_comp),2.)+(2./7.)*pow(abs(D1pRe_comp),2.)-(4./7.)*pow(abs(D2pRe_comp),2.),0.);
        pair <double,double> H1_21((2./sqrt(5.))*real(S0pRe_comp*conj(D1pRe_comp))+(2.*sqrt(3.)/5.)*real(P0pRe_comp*conj(P1pRe_comp))+(2./7.)*real(D0pRe_comp*conj(D1pRe_comp)),0.);
	pair <double,double> H0_21(H1_21.first+(2./7.*sqrt(6.))*real(D1pRe_comp*conj(D2pRe_comp)),0.);
        pair <double,double> H0_22((2./sqrt(5.))*real(S0pRe_comp*conj(D2pRe_comp))-(4./7.)*real(D0pRe_comp*conj(D2pRe_comp)),0.);
        pair <double,double> H1_22(H0_22.first+(sqrt(6.)/7.)*pow(abs(D1pRe_comp),2.)+(sqrt(6.)/5.)*pow(abs(P1pRe_comp),2.),0.);




	outfile << H1_00.first << "\t"<< H1_00.second <<"\t" << H0_00.first <<"\t" << H0_00.second  <<"\t"
		<< H1_10.first << "\t"<< H1_10.second <<"\t" << H0_10.first <<"\t" << H0_10.second  <<"\t"
		<< H1_11.first << "\t"<< H1_11.second <<"\t" << H0_11.first <<"\t" << H0_11.second  <<"\t"
		<< H1_20.first << "\t"<< H1_20.second <<"\t" << H0_20.first <<"\t" << H0_20.second  <<"\t"
		<< H1_21.first << "\t"<< H1_21.second <<"\t" << H0_21.first <<"\t" << H0_21.second  <<"\t"
		<< H0_22.first << "\t"<< H0_22.second <<"\t" << H1_22.first <<"\t" << H1_22.second  <<"\t";

        outfile << endl;

	//	if( i==44 && j==3) cout <<" H0_11= "<<H0_11.first<<"   "<<H1_11.first<<endl;

        
        chdir( ".." );
    }
    }
    return 0;
}



