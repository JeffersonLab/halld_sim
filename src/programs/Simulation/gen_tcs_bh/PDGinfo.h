#ifndef PDGinfo_h
#define PDGinfo_h
#include <iostream>
#include <Riostream.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
using namespace std;

/*
// particle names following GEANT particle numbering
string GeantPID_name[49]={"0","gamma","positron","electron","neutrino","mu+","mu-","pi0","pi+","pi-",
	"K0long","K+","K-","neutron","proton","antiproton","K0short","eta","lambda","sigma+",
	"sigma0","sigma-","xi0","xi-","omega","antineutron","antilambda","antisigma-","antisigma0","antisigma+",
	"antixi0","antixi+","antiOmega+","tau+","tau-","D+","D-","D0","antiD0","f+",
	"f-","lambdac+","W+","W-","Z0","deuteron","tritium","alpha","genatino"};
*/
const double GeantPID_mass[49]={0,0,0.00051099891,0.00051099891,0,0.105658367,0.105658367,0.1349766,0.13957018,0.13957018,
	0.497614,0.493677,0.493677,0.93956546,0.93827203,0.93827203, 0.497614,0.547853,1.115683,1.18937,
	1.192642,1.197449,1.31486,1.32171,1.67245,0.93956546,1.115683,1.197449,1.192642,1.18937,
	1.31486,1.32171,1.67245,1.77684,1.77684,1.86962,1.86962,1.86484,1.86484,
	0,0,0,0,0,0,0,0,0,0};
	//,"f+","f-","lambdac+","W+","W-","Z0","deuteron","tritium","alpha","genatino");

// X0 in g/cm2 calculated ordered by element
//double PDG_radlenght_Z[93]={0,63.04,125.97,94.32,82.78,  };
const double RadLenght[63]=
	{	30050,
		890.4, 125.97/0.169, 94.32/0.125, 94.32/0.125,79.61,19.32,40.87, 82.78/0.534 ,65.19/1.848,
		28.93/1.204,30.390, 42.7/3.52, 8.897, 21.82/2.329 , 37.99/0.807, 34.24/1.141, 19.28/1.574, 32.93/1.507,2.05,
		36.1, 2.59,16.16/4.54, 1.265,0.5612, 19.55/1.396, 1.757, 26.54/2.65 , 6./18.95 ,1.43,
		42.9, 42.4, 12.25/5.323, 8.82/7.31, 8.48/2.953, 1.12, 30.39, 6.76/19.3, 6.54/21.45, 6.46/19.32,
		27.94/3.97,36.2/1.842, 36.2/1.563, 1.243/4.51, 39.26/2.635, 7.39/8.3, 27.05/2.2, 21.91/2.17, 26.57/2.3, 28.17/2.23,
		46.47/0.667, 46.56/1.263, 45.23/2.489, 45./0.703,44.85/0.93, 41.92/1.18, 41.5/1.2, 44.77/0.89, 39.95/1.4, 40.55/1.19,
		44.77/0.9, 1.671/2.2, 1.956/1.03	
	};

class PDGinfo{
	public:
		string TargetMat(int);
		string PDG_PID_name(int);
		string GeantPID_name(int);
	private:
		vector <string> fGeantPID_name();
		vector <string> fTargetMat();
		vector <string> fPDG_PID_name();
}; 

#endif
