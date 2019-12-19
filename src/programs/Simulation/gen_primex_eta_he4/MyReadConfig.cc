/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2019       GlueX and PrimEX-D Collaborations               * 
*                                                                         *                                                                                                                               
* Author: The GlueX and PrimEX-D Collaborations                           *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#include "MyReadConfig.h"

//------------------------------------------------------------------------------
MyReadConfig::MyReadConfig() {
  nLine=0;
  
  // Check if env. variable $ exist
  //
  
  strCaLibPath = gSystem->pwd();
  
  this->ReadConfigFile( "config.dat" );
}

//------------------------------------------------------------------------------
MyReadConfig::MyReadConfig( Char_t * szFin ) {
  nLine=0;
  
  this->ReadConfigFile( szFin );
}

//------------------------------------------------------------------------------
MyReadConfig::~MyReadConfig() {
}

//------------------------------------------------------------------------------
void MyReadConfig::ReadConfigFile( const Char_t *szFin ) {
  // Build File name
  char szCalibFile[128];
  sprintf( szCalibFile, 
	   "%s/%s",
	   strCaLibPath.Data(),
	   szFin );
  
  //
  ifstream infile;
  infile.open( szCalibFile );
  
  if ( !infile.is_open() ) {
      printf("\n ERROR: opening \"%s\" file ! ! !\n\n", szFin);
  } else {
      printf("\n ---------------------------------------- \n");
      printf(" Read File : \"%s\"\n", szFin);
      
      while ( infile.good() ) {
	strLine[nLine]="";
	strLine[nLine].ReadLine(infile);
	
	if ( strLine[nLine].BeginsWith("#") ) {
	  // 	      printf("Comment : %s \n", strLine[nLine].Data());
	} else if ( strLine[nLine].Contains("FILE") ) {
	  TString strFile = this->ExtractName( strLine[nLine] );
	  
	  this->ReadConfigFile( strFile.Data() );
	}
	nLine++;
	if ( nLine >= MAX_LINE ) {
	  printf("\n ERROR: Number of lines is more than MAX_LINE ! ! !\n\n");
	  gSystem->Exit(0);
	}
      }
  }
  
  infile.close();
  return;
}

//------------------------------------------------------------------------------
TString MyReadConfig::ExtractName( TString strIn ) {
  Ssiz_t aa = strIn.First(":")+1;
  Ssiz_t bb = strIn.Length()-aa;
  
  TString cc = strIn(aa, bb);
  cc.ReplaceAll(" ","");
  return cc;
}

//------------------------------------------------------------------------------
TString MyReadConfig::GetConfigName( TString name ) {
  TString strOut;
  
  for ( Int_t i=0; i<nLine; i++) {
    if ( !( strLine[i].BeginsWith("#") ) && 
	  strLine[i].Contains( name )    ) {
      strOut = this->ExtractName( strLine[i] );
    }      
  }
  return strOut; 
}

//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig6Par( TString name ) {
  Char_t tmp1[256];
  Char_t tmp2[256];
  Char_t tmp3[256];
  Char_t tmp4[256];
  Char_t tmp5[256];
  Char_t tmp6[256];
  
  Double_t* iPar = new Double_t[6];

  for ( Int_t i=0; i<nLine; i++) {
    if ( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
      sscanf( strLine[i].Data(),"%*s %s %s %s %s %s %s",tmp1,tmp2,tmp3,tmp4,tmp5,tmp6);
    }     
  }
  
  iPar[0] = atof(tmp1);
  iPar[1] = atof(tmp2);
  iPar[2] = atof(tmp3);
  iPar[3] = atof(tmp4);
  iPar[4] = atof(tmp5);
  iPar[5] = atof(tmp6);
  
  return iPar; 
}

//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig5Par( TString name ) {
  Char_t tmp1[256];
  Char_t tmp2[256];
  Char_t tmp3[256];
  Char_t tmp4[256];
  Char_t tmp5[256];
  
  Double_t* iPar = new Double_t[5];

  for ( Int_t i=0; i<nLine; i++) {
    if ( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
      sscanf( strLine[i].Data(),"%*s %s %s %s %s %s",tmp1,tmp2,tmp3,tmp4,tmp5);
    }     
  }
  
  iPar[0] = atof(tmp1);
  iPar[1] = atof(tmp2);
  iPar[2] = atof(tmp3);
  iPar[3] = atof(tmp4);
  iPar[4] = atof(tmp5);
  
  return iPar; 
}

//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig4Par( TString name ) {
  Char_t tmp1[256];
  Char_t tmp2[256];
  Char_t tmp3[256];
  Char_t tmp4[256];
  
  Double_t* iPar = new Double_t[4];

  for ( Int_t i=0; i<nLine; i++) {
    if ( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
      sscanf( strLine[i].Data(),"%*s %s %s %s %s",tmp1,tmp2,tmp3,tmp4);
    }     
  }
  
  iPar[0] = atof(tmp1);
  iPar[1] = atof(tmp2);
  iPar[2] = atof(tmp3);
  iPar[3] = atof(tmp4);
  
  return iPar; 
}

//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig3Par( TString name ) {
  Char_t tmp1[256];
  Char_t tmp2[256];
  Char_t tmp3[256];
  
  Double_t* iPar = new Double_t[3];

  for ( Int_t i=0; i<nLine; i++) {
      if ( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
	sscanf( strLine[i].Data(),"%*s %s %s %s",tmp1,tmp2,tmp3);
      }     
  }
  
  iPar[0] = atof(tmp1);
  iPar[1] = atof(tmp2);
  iPar[2] = atof(tmp3);
  
  return iPar; 
}

//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig2Par( TString name ) {
  Char_t tmp1[256];
  Char_t tmp2[256];
  
  Double_t* iPar = new Double_t[2];

  for ( Int_t i=0; i<nLine; i++) {
    if( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
      sscanf( strLine[i].Data(),"%*s %s %s",tmp1,tmp2);
    }     
  }
  
  iPar[0] = atof(tmp1);
  iPar[1] = atof(tmp2);
   
  return iPar; 
}




//------------------------------------------------------------------------------

Double_t* MyReadConfig::GetConfig1Par( TString name ) {
  Char_t tmp1[256];
  
  Double_t* iPar = new Double_t[1];
  
  for ( Int_t i=0; i<nLine; i++) {
    if ( !( strLine[i].BeginsWith("#") ) && strLine[i].Contains( name )    ) {
      sscanf( strLine[i].Data(),"%*s %s",tmp1);
    }     
  }
  
  iPar[0] = atof(tmp1);
  
  return iPar; 
}


