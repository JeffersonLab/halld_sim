#ifndef FillTables_h
#define FillTables_h
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "TableProvider.h"
#include "TabUtils.h"
#include "Options_tcs.h"
#include "TreatOptions.h"
#include "Utils.h"
#include <unistd.h>
#include "Constants.h"
#include "genSettings_t.h"

class FillTables {

public:
int FillTable(int reaction, TableProvider& table_provider);

private:

};

#endif



