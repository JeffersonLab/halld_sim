#define PDGinfo_cxx
#include "PDGinfo.h"
using namespace std;
extern const double GeantPID_mass[49],  RadLenght[63];
//extern string GeantPID_name[49], PDG_PID_name[560], TargetMat[63];

vector <string> PDGinfo::fTargetMat(){
	ifstream infile; string line;
	string ff = "../Param/Materials.txt";
	vector <string> TM;
	infile.open("ff");
	if (!infile){
		TM.push_back( "ErrornoMat");
	} 
	while(infile.good()){
		getline(infile,line);
		TM.push_back(line);
	}	
	return TM;

}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

vector <string> PDGinfo::fGeantPID_name(){
	
	ifstream infile; string line;
	string ff = "../Param/Geant_names.txt";
	vector <string> TM;
	infile.open("ff");
	if (!infile){
		TM.push_back( "Errornoname");
	} 
	while(infile.good()){
		getline(infile,line);
		TM.push_back(line);
	}
	
	return TM;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

vector <string> PDGinfo::fPDG_PID_name(){
	
	ifstream infile; string line;
	string ff = "../Param/PDGnames.txt";
	vector <string> TM;
	infile.open("ff");
	if (!infile){
		TM.push_back( "ErrornoMat");
	} 
	while(infile.good()){
		getline(infile,line);
		TM.push_back(line);
	}
	for (int i=0;i<2250;i++){
		if (i==990) TM.push_back("pomeron");
		else if (i==1114) TM.push_back("Delta-");
		else if (i==2112) TM.push_back("neutron");
		else if (i==2114) TM.push_back("Delta0");
		else if (i==2212) TM.push_back("proton");
		else if (i==2214) TM.push_back("Delta+");
		else if (i==2224) TM.push_back("Delta++");
		else TM.push_back(to_string(i));
	}
	
	return TM;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

string PDGinfo::TargetMat(int i){
	vector <string> TM=fTargetMat();
  	if (TM.size()>(unsigned int) i) return TM[i];
	else return "noMat";
} 

//-------------------------------------------------------------------
//-------------------------------------------------------------------

string PDGinfo::PDG_PID_name(int i){
	vector <string> TM=fPDG_PID_name();
  	if (TM.size()>(unsigned int) i) return TM[i];
	else return "noname";
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

string PDGinfo::GeantPID_name(int i){
	vector <string> TM=fGeantPID_name();
	if (TM.size()>(unsigned int) i) return TM[i];
        else return "noname";
}
 
//-------------------------------------------------------------------
//-------------------------------------------------------------------
