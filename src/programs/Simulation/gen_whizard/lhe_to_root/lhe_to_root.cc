/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include "Riostream.h"
#include <TH3.h>
#include <TH2.h>

int main(int argc, char * argv[]) {
  
  TString RootDirName = "";
  if (argc == 2) {
    RootDirName = std::string(argv[1]);
  } else {
    std::cout << "lhe_to_root" << std::endl;
    std::cout << "===========" << std::endl;
    std::cout << "This program converts a LHE file into a ROOT file containing the events in the lhe tree" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage:                                                                                 " << std::endl;
    std::cout << "Method : ./lhe_to_root input_directory_path/                                           " << std::endl;
  }
  
  TSystemDirectory dir(RootDirName, RootDirName);
  TList * files = dir.GetListOfFiles(); 
 
  Char_t tmp0[256];
  Char_t tmp1[256];
  Char_t tmp2[256];
  Char_t tmp3[256];
  Char_t tmp4[256];
  Char_t tmp5[256];
  Char_t tmp6[256];
  Char_t tmp7[256];
  Char_t tmp8[256];
  
  int npart = 0;
  int pdg[100];
  double xs = 0, er_xs = 0, weight = 0;
  double px[100], py[100], pz[100], e[100], m[100];
  int status[100], first_daughter[100], last_daughter[100];
  for (int i = 0; i < 100; i ++) {
    status[i] = 0;
    first_daughter[i] = 0;
    last_daughter[i] = 0;
    pdg[i] = 0;
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
    e[i] = 0;
    m[i] = 0;
  }
  
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next()) ) {
      fname = file->GetName();
      //if (!file->IsDirectory() && fname.Contains("gz")) {
      if (!file->IsDirectory() && fname.Contains("lhe")) {
	//cout << fname << endl;
	TString GZipFileName = RootDirName + fname;
	TString FileUnZip = GZipFileName;
	//FileUnZip.ReplaceAll(".gz","");
	//system(TString::Format("gunzip %s", GZipFileName.Data()));
	std::cout << GZipFileName << std::endl;
	TString RootFileName = GZipFileName;
	//RootFileName.ReplaceAll(".lhe.gz",".root");
	RootFileName.ReplaceAll(".lhe",".root");
	TFile * m_file = new TFile(RootFileName, "RECREATE");
	TTree * m_tree = new TTree("lhe", "LHE tree w/ all particles, weight, and xs [mb]");
	m_tree->Branch("npart",&npart,"npart/I");
	m_tree->Branch("weight",&weight,"weight/D");
	m_tree->Branch("xs",&xs,"xs/D");
	m_tree->Branch("er_xs",&er_xs,"er_xs/D");
	m_tree->Branch("pdg", pdg, "pdg[npart]/I");
	m_tree->Branch("status", status, "status[npart]/I");
	m_tree->Branch("first_daughter", first_daughter, "first_daughter[npart]/I");
	m_tree->Branch("last_daughter", last_daughter, "last_daughter[npart]/I");
	m_tree->Branch("px", px, "px[npart]/D");
	m_tree->Branch("py", py, "py[npart]/D");
	m_tree->Branch("pz", pz, "pz[npart]/D");
	m_tree->Branch("e", e, "e[npart]/D");
	m_tree->Branch("m", m, "m[npart]/D");
	
	ifstream in;
	int l = 0;
	int k = 0;
	in.open(FileUnZip);
	int j = 0;
	while(in.good()) {
	  TString LogLine = "";
	  LogLine.ReadLine(in);
	  j ++;
	  if (j == 6) {
	    sscanf(LogLine.Data(),"%s %s %*s %*s", tmp0, tmp1);
	    xs = atof(tmp0); //pb
	    xs *= 1e-9;//pb to mb
	    er_xs = atof(tmp1); //pb
	    er_xs *= 1e-9;//pb to mb
	  }
	  /*
	  if(LogLine.Contains("xsecinfo neve")) {
	    //cout << LogLine << endl;
	    LogLine.ReplaceAll("\"","");
	    LogLine.ReplaceAll("="," ");
	    sscanf(LogLine.Data(),"%*s %*s %*s %*s %s",tmp0);
	    //cout << "char nevt= " << tmp0 << " xs= " << tmp1 << endl; 
	    //cout << "numb nevt= " << atof(tmp0) << " xs= " << atof(tmp1) << endl; 
	    xs = atof(tmp0); //pb
	    xs *= 1e-9;//pb to mb
	    //cout <<"xs " << xs << " mb " << endl;
	  }
	  */
	  if (LogLine.Contains("weight name")) {
	    LogLine.ReplaceAll("\"","");
            LogLine.ReplaceAll("="," ");
	    LogLine.ReplaceAll(">"," ");
	    LogLine.ReplaceAll("<"," ");
	    //cout << LogLine << endl; 
	    sscanf(LogLine.Data(),"%*s %*s %*s %s %*s",tmp0);
	    weight = atof(tmp0);
	    //cout <<"char weight = " << tmp0 << endl; 
	    //cout <<"numb weight = " << atof(tmp0) << endl; 
	    //cout <<"weight " << weight << endl;
	  }
	  //in>>npart>>dummy>>dummy>>dummy>>dummy>>dummy;
	  if (l == 1) {
	    //cout << LogLine << endl;
	    sscanf(LogLine.Data(),"%s %*s %*s %*s %*s %*s",tmp0);
	    npart = atof(tmp0);
	    k = 0;
	  }
	  if (l < (npart + 2) && l > 1) {
	    //11 1 1 2 0 0  2.7259169828E-02 -1.9383361874E-02  8.7489276154E+00  8.7489915681E+00  5.1100000000E-04  0.0000000000E+00  9.0000000000E+00
	    sscanf(LogLine.Data(),"%s %s %s %s %*s %*s %s %s %s %s %s %*s %*s", tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8);
	    pdg[k] = atof(tmp0);
	    status[k] = atof(tmp1);
	    first_daughter[k] = atof(tmp2);
	    last_daughter[k] = atof(tmp3);
	    px[k] = atof(tmp4);
	    py[k] = atof(tmp5);
	    pz[k] = atof(tmp6);
	    e[k] = atof(tmp7);
	    m[k] = atof(tmp8);
	    k++;
	    //cout << LogLine << endl;
	  }
	  if (LogLine.Contains("<event>")) {
	    l = 0;
	  }
	  l ++;
	  if (LogLine.Contains("</event>")) {
	    //npart = 0; 
	    //cout << "npart " << npart << endl;
	    //cout <<"l = " << l << endl;
	    //l = 0;
	    m_tree->Fill();
	    //cout <<k<<endl;
	  }
	}
	m_file->Write();
	m_file->Close();
	delete m_file;
	m_tree = NULL;
	//system(TString::Format("gzip %s", FileUnZip.Data()));
      }
    }
  }
}
