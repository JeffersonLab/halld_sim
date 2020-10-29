/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>
#include <complex>
#include <string>
#include <time.h>

#include "IUAmpTools/FitResults.h"
#include "TFile.h"

#include "wave.h"
#include "moment.h"

#include "TFile.h"





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
    


   //! Set waveset, has to be same order as in fit.cfg!

 // vector<wave> negative;
  //negative.push_back(wave("S0", 0, 0));

  vector<wave> positive;
  positive.push_back(wave("S0+", 0, 0));
  positive.push_back(wave("D0+", 2, 0));
  positive.push_back(wave("D1+", 2, 1));
 positive.push_back(wave("D2+", 2, 2));
  positive.push_back(wave("P0+", 1, 0));
   positive.push_back(wave("P1+", 1, 1));
positive.push_back(wave("G0+", 4, 0));
positive.push_back(wave("G1+", 4, 1));

  coherent_waves wsPos, wsNeg;
  wsPos.reflectivity = +1;
  wsPos.waves = positive;

  //wsNeg.reflectivity = -1;
  //wsNeg.waves = negative;

  waveset ws;
  ws.push_back(wsPos);
  //ws.push_back(wsNeg);
  
  //take for index step size 4 as there are two of the same waves next to each other corresponding to diff. sums
  size_t lastIdx = 0;
  for (size_t i = 0; i < ws.size(); i++) //ws.size gives number of coherent sums (=2 in this case, negative, positive)
    for (size_t j = 0; j < ws[i].waves.size(); j++, lastIdx += 4)// ws[i].waves.size() gives number of waves in given sum, index is increased by two for next wave as each wave takes two index for real and imaginary components
      ws[i].waves[j].setIndex(lastIdx);
    

  //LMAX=2*l_max (l_max is highest wave)
  //We consider M values form 0 to L, as H(LM)=H(L-M) due to parity invariance
  size_t LMAX;    //highest wave
  Biggest_lm(ws, &LMAX);


  double step = ( highMass - lowMass ) / kNumBins;
    double stept = ( hight - lowt ) / kNumBinst;
    
    ofstream outfile;
    outfile.open( outfileName.c_str() );

    outfile <<"M"<<"\t"<<"t"; //First line contains names of variables, first two colomns correspond to M(invariant mass) and t

    for (int L = 0; L<= int(LMAX); L++) {// the rest of the colomn correspond to moments
      for (int M = 0; M<= L; M++) {
    
	outfile<<"\t"<<"H0_"<<L<<M<<"\t"<<"H0_"<<L<<M<<"uncert.";
	outfile<<"\t"<<"H1_"<<L<<M<<"\t"<<"H1_"<<L<<M<<"uncert.";
      }}
    outfile<<endl;



    // descend into the directory that contains the bins
    chdir( fitDir.c_str() );



    //Looping through M_eta_pi and t bins
    for( int i = 0; i < kNumBins;i++ ){
      for( int j = 0; j < kNumBinst; j++ ){  
	cout<<"bin "<<i<<"_"<<j<<endl;
        
        ostringstream dir;
        dir << "bin_" << i<<"_"<<j;
        chdir( dir.str().c_str() );
        
        string resultsFile;
        resultsFile ="bin_" + std::to_string(i)+"_"+ std::to_string(j)+".fit";
        FitResults results(resultsFile);
	

	if( !results.valid() ){
	  outfile<< lowMass + step * i + step / 2. << "\t";
	  outfile<< lowt + stept * j + stept / 2. << "\t";
	  for (int L = 0; L<= pow(LMAX,1); L++) {// calculating moments and writing to a file
	    for (int M = 0; M<= L; M++) {
     
	      outfile << 0<< "\t"<< 0 <<"\t";
	      outfile << 0<< "\t"<< 0 <<"\t";
	    }}
	  outfile << endl;

	  chdir( ".." );   
	  continue;
        }



        
        // print out the bin center
        outfile << lowMass + step * i + step / 2. << "\t";
	outfile << lowt + stept * j + stept / 2. << "\t";
  
  
    for (int L = 0; L<= pow(LMAX,1); L++) {// calculating moments and writing to a file
    for (int M = 0; M<= L; M++) {
     
      outfile << real(decomposeMoment(0, L, M, ws, results.parValueList()))<< "\t"<< 0 <<"\t";
      outfile << real(decomposeMoment(1, L, M, ws, results.parValueList()))<< "\t"<< 0 <<"\t";
      // if(L==0 && M==0 && i==44 && j==3)cout <<" H0_00=  "<< real(decomposeMoment(0, L, M, ws, results.parValueList()))<<" H1_00=  "<< real(decomposeMoment(1, L, M, ws, results.parValueList()))<<" H2_00=  "<< decomposeMoment(2, L, M, ws, results.parValueList())<<endl;

        }}
        outfile << endl;
        
        chdir( ".." );
    }
    }




    return 0;
}
