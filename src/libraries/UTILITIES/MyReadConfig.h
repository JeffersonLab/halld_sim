/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2019       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/


#if !defined(MYREADCONFIG)
#define MYREADCONFIG

//#ifndef _MYREADCONFIG_H_
//#define _MYREADCONFIG_H_


#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include <string>
#include <math.h>
#include "TFile.h"
#include "TObject.h"
#include "TSystem.h"
#include "TLorentzVector.h"

const Int_t MAX_LINE = 1000;

class MyReadConfig {
 private:
  Int_t   nLine;
  TString strLine[MAX_LINE];

  TString strCaLibPath;

  //void ReadConfigFile( const Char_t* );
  
 protected:
 public:
  MyReadConfig();
  MyReadConfig( Char_t * );
  virtual ~MyReadConfig();

  TString ExtractName(TString);
  
  TString GetConfigName( TString );
  void ReadConfigFile( const Char_t* );
  

  Double_t* GetConfig1Par( TString );
  Double_t* GetConfig2Par( TString );
  Double_t* GetConfig3Par( TString );
  Double_t* GetConfig4Par( TString );
  Double_t* GetConfig5Par( TString );
  Double_t* GetConfig6Par( TString );

};

#endif // _MYREADCONFIG_H_
