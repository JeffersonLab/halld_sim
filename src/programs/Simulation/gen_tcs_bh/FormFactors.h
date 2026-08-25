#ifndef FormFactors_h
#define FormFactors_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "Options_tcs.h"
#include "TreatOptions.h"
#include "Utils.h"
#include <unistd.h>
#include "Constants.h"

float G_E(int,float);
float G_M(int,float);
float FormFactors(int FF,int nucleon, float QQ);

#endif



