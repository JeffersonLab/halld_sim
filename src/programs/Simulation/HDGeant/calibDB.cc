
#include <stdlib.h>
#include <dlfcn.h>
#include <unistd.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include <JANA/Calibrations/JCalibrationFile.h>
#include <JANA/Calibrations/JCalibrationManager.h>
#include <DANA/DApplication.h>
#include <HDGEOMETRY/DGeometry.h>
#include <HDGEOMETRY/DMagneticFieldMapFineMesh.h>
#include <HDGEOMETRY/DMagneticFieldMapCalibDB.h>
#include <HDGEOMETRY/DMagneticFieldMapConst.h>
#include "HDGEOMETRY/DMagneticFieldMapSpoiled.h"
#include "HDGEOMETRY/DMagneticFieldMapParameterized.h"
#include "HDGEOMETRY/DMagneticFieldMapNoField.h"
#include "HDGEOMETRY/DMagneticFieldMapPSConst.h"
#include "HDGEOMETRY/DMagneticFieldMapPS2DMap.h"

extern "C" {
#include "calibDB.h"
};
#include "controlparams.h"


extern "C" int hddsgeant3_runtime_(void);  // called from uginit.F. defined in calibDB.cc
extern "C" void md5geom_(char *md5);
void init_runtime_xml(void);
void md5geom_runtime(char *md5);
extern "C" const char* GetMD5Geom(void);

extern "C" {
   int getcalib_(const char* namepath, unsigned int *Nvals, float* vals) {
      int retval;
      char name[999];
      int n;
      for (n=0; n < 999 && namepath[n] != ' '; ++n)
         name[n] = namepath[n];
      name[n] = 0;
      retval = GetCalib(name, Nvals, vals);
      return retval;
   }
};

bool nofield=false;
DMagneticFieldMap *Bmap=NULL;
DMagneticFieldMapPS *PS_Bmap=NULL;
static JCalibration *jcalib=NULL;

extern "C" {
   void md5geom_wrapper_(char *md5);
}


//----------------
// initcalibdb_
//----------------
void initcalibdb_(char *bfield_type, char *bfield_map, char *PS_bfield_type, char *PS_bfield_map, int *runno)
{
   ios::sync_with_stdio(true);

   if(!japp){
     _DBG_<<" JApplication missing, exiting !!"<<endl;
     exit(-1) ;
   }

   // Get the JCalibration object
   jcalib =  japp->GetService<JCalibrationManager>()->GetJCalibration(*runno);
 
   // The actual DMagneticFieldMap subclass can be specified in
   // the control.in file. Since it is read in as integers of
   // "MIXED" format through ffkey though (who knows what that
   // means!) then there can be trailing white space at the end
   // of the string. Here, we terminate the string with a null
   // to eliminate that white space.
   while(strlen(bfield_type)>0 && bfield_type[strlen(bfield_type)-1]==' ')bfield_type[strlen(bfield_type)-1] = 0;
   while(strlen(bfield_map)>0 && bfield_map[strlen(bfield_map)-1]==' ')bfield_map[strlen(bfield_map)-1] = 0;
   while(strlen(PS_bfield_type)>0 && PS_bfield_type[strlen(PS_bfield_type)-1]==' ')PS_bfield_type[strlen(PS_bfield_type)-1] = 0;
   while(strlen(PS_bfield_map)>0 && PS_bfield_map[strlen(PS_bfield_map)-1]==' ')PS_bfield_map[strlen(PS_bfield_map)-1] = 0;
   
   // Read in the field map from the appropriate source
   if(bfield_type[0] == 0)strcpy(bfield_type, "CalibDB");
   string bfield_type_str(bfield_type);

   const char *ccdb_help =
	   " \n"
	   " Could not load the solenoid field map from the CCDB!\n"
	   " Please specify the solenoid field map to use on the command line, e.g.:\n"
	   " \n"
	   "   -PBFIELD_MAP=Magnets/Solenoid/solenoid_1200A_poisson_20140520\n"
	   " or\n"
	   "   -PBFIELD_TYPE=NoField\n";

   if(bfield_type_str=="CalibDB"){
     // if the magnetic field is specified in control.in, then use that value instead of the CCDB values
     if(strlen(bfield_map))
       Bmap = new DMagneticFieldMapFineMesh(japp,*runno,bfield_map);
     else {
          
       // see if we can load the name of the magnetic field map to use from the calib DB
       map<string,string> bfield_map_name;
       if(jcalib->GetCalib("/Magnets/Solenoid/solenoid_map", bfield_map_name)) {
	 // if we can't find information in the CCDB, then quit with an error message
	 _DBG_<<ccdb_help<<endl;
	 exit(-1);
       } else {
	 if( bfield_map_name.find("map_name") != bfield_map_name.end() ) {
	   if( bfield_map_name["map_name"] == "NoField" )     // special case for no magnetic field
	     Bmap = new DMagneticFieldMapNoField(japp);
	   else
	     Bmap = new DMagneticFieldMapFineMesh(japp,*runno,bfield_map_name["map_name"]);  // pass along the name of the magnetic field map to load
	 } else {
	   // if we can't find information in the CCDB, then quit with an error message
	   jerr << ccdb_help << endl;
	   exit(-1);
	 }
       }
     }
   }
   else if(bfield_type_str=="NoField"){
     nofield=true;
   }
   else if(bfield_type_str=="Const"){
     Bmap = new DMagneticFieldMapConst(0.0, 0.0, 1.9);
   }
   else{
     _DBG_<<" Unknown DMagneticFieldMap subclass \"DMagneticFieldMap"<<bfield_type_str<<"\" !!"<<endl;
     exit(-1);
   }   

   // also load the PS magnet field map similarly to the solenoid field map above
   if(PS_bfield_type[0] == 0)strcpy(PS_bfield_type, "CalibDB");
   string PS_bfield_type_str(PS_bfield_type);

   const char *PS_ccdb_help =
	   " \n"
	   " Could not load the pair spectrometer field map from the CCDB!\n"
	   " Please specify the pair spectrometer field map to use on the command line, e.g.:\n"
	   " \n"
	   "   -PBFIELD_MAP=Magnets/PairSpectrometer/PS_1.8T_20150513_test\n"
	   " or\n"
	   "   -PBFIELD_TYPE=Const\n";

   if(PS_bfield_type_str=="CalibDB"){
     // if the magnetic field is specified in control.in, then use that value instead of the CCDB values
     if(strlen(PS_bfield_map))
       PS_Bmap = new DMagneticFieldMapPS2DMap(japp,*runno,PS_bfield_map);
     else {
          
       // see if we can load the name of the magnetic field map to use from the calib DB
       map<string,string> PS_bfield_map_name;
       if(jcalib->GetCalib("/Magnets/PairSpectrometer/ps_magnet_map", PS_bfield_map_name)) {
	 // if we can't find information in the CCDB, then quit with an error message
	 _DBG_<<PS_ccdb_help<<endl;
	 exit(-1);
       } else {
	 if( PS_bfield_map_name.find("map_name") != PS_bfield_map_name.end() ) {
	   PS_Bmap = new DMagneticFieldMapPS2DMap(japp,*runno,PS_bfield_map_name["map_name"]);  // pass along the name of the magnetic field map to load
	 } else {
	   // if we can't find information in the CCDB, then quit with an error message
	   jerr << PS_ccdb_help << endl;
	   exit(-1);
	 }
       }
     }
   }
   else if(PS_bfield_type_str=="Const"){
     PS_Bmap = new DMagneticFieldMapPSConst(0.0, 1.64, 0.0);
   }
   else{
     _DBG_<<" Unknown DMagneticFieldMap subclass \"DMagneticFieldMap"<<PS_bfield_type_str<<"\" !!"<<endl;
     exit(-1);
   }   

}

//----------------
// gufld_ps_
//----------------
void gufld_ps_(float *r, float *B)
{
   /// Wrapper function to allow the FORTRAN gufld routine to
   /// use the C++ class DMagneticFieldMap to access the 
   /// B-field.

   if(!PS_Bmap){
      _DBG_<<"Call to gufld_ps when PS_Bmap not intialized! Exiting."<<endl;
      exit(-1);
   }
   
   double x = r[0];
   double y = r[1];
   double z = r[2];
   double Bx, By, Bz;
   
   PS_Bmap->GetField(x, y, z, Bx, By, Bz);

   B[0] = Bx;
   B[1] = By;
   B[2] = Bz;
}

//----------------
// gufld_db_
//----------------
void gufld_db_(float *r, float *B)
{
   /// Wrapper function to allow the FORTRAN gufld routine to
   /// use the C++ class DMagneticFieldMap to access the 
   /// B-field.
  if (nofield){
    B[0]=0.;
    B[1]=0.;
    B[2]=0.;
    
    return;
  }
  

   if(!Bmap){
      _DBG_<<"Call to gufld_db when Bmap not intialized! Exiting."<<endl;
      exit(-1);
   }
   
   double x = r[0];
   double y = r[1];
   double z = r[2];
   double Bx, By, Bz;
   
   Bmap->GetField(x, y, z, Bx, By, Bz);

   B[0] = Bx;
   B[1] = By;
   B[2] = Bz;
}

//----------------
// GetCalib
//----------------
int GetCalib(const char* namepath, unsigned int *Nvals, float* vals)
{
   /// C-callable routine for accessing calibration constants.
   /// The values specified by "namepath" will be read into the array
   /// "vals". The "vals" array should have enough memory allocated
   /// to hold *Nvals elements. If not, only the first *Nvals elements
   /// will be copied and a non-zero value returned. If the number
   /// of values in the database are less than *Nvals, then all values
   /// are copied, *Nvals is updated to reflect the number of valid
   /// elements in "vals", and a value of 0 is returned.
   
   if(!jcalib){
      _DBG_<<"ERROR - GetCalib() called when jcalib not set!"<<endl;
      _DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }

   vector<float> vvals;
   jcalib->Get(namepath, vvals);
   if(vvals.size()<*Nvals)*Nvals = vvals.size();
   for(unsigned int i=0; i<*Nvals; i++)vals[i] = vvals[i];
   
   return vvals.size()>*Nvals; // return 0 if OK, 1 if not
}

//----------------
// GetLorentzDeflections
//----------------
void GetLorentzDeflections(float *lorentz_x, float *lorentz_z, float **lorentz_nx, float **lorentz_nz
   , const unsigned int Nxpoints, const unsigned int Nzpoints)
{
   /// C-callable routine for accessing calibration constants.
   /// The values specified by "namepath" will be read into the array
   /// "vals". The "vals" array should have enough memory allocated
   /// to hold *Nvals elements. If not, only the first *Nvals elements
   /// will be copied and a non-zero value returned. If the number
   /// of values in the database are less than *Nvals, then all values
   /// are copied, *Nvals is updated to reflect the number of valid
   /// elements in "vals", and a value of 0 is returned.
   
   // Make sure jcalib is set
   if(!jcalib){
      _DBG_<<"ERROR - GetLorentzDeflections() called when jcalib not set!"<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }

   // Get constants and do basic check on number of elements
   vector< map<string, float> > tvals;
   jcalib->Get("FDC/lorentz_deflections", tvals);
   if(tvals.size() != Nxpoints*Nzpoints){
      _DBG_<<"ERROR - GetLorentzDeflections() number of elements in calib DB"<<endl;
      _DBG_<<"ERROR - not the same as expected. DB="<<tvals.size()<<" expected"<<Nxpoints*Nzpoints<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }
   
   // Notify user
   cout<<"Read "<<tvals.size()<<" values from FDC/lorentz_deflections in calibDB"<<endl;
   cout<<"   lorentz_deflections columns (alphabetical): ";
   map<string,float>::iterator iter;
   for(iter=tvals[0].begin(); iter!=tvals[0].end(); iter++)cout<<iter->first<<" ";
   cout<<endl;
   
   // Copy values into tables. We preserve the order since that is how it was
   // originally done in hitFDC.c
   for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      unsigned int xindex = i/Nzpoints;
      unsigned int zindex = i%Nzpoints;
      lorentz_x[xindex] = row["x"];
      lorentz_z[zindex] = row["z"];
      lorentz_nx[xindex][zindex] = row["nx"];
      lorentz_nz[xindex][zindex] = row["nz"];
   }   
}
//----------------
// GetConstants
//----------------
int GetConstants(const char* namepath, int *Nvals, float* vals, mystr_t *strings)
{
   /// C-callable routine for accessing calibration constants.
   /// The values specified by "namepath" will be read into the array
   /// "vals". The "vals" array should have enough memory allocated
   /// to hold *Nvals elements. If not, only the first *Nvals elements
   /// will be copied and a non-zero value returned. If the number
   /// of values in the database are less than *Nvals, then all values
   /// are copied, *Nvals is updated to reflect the number of valid
   /// elements in "vals", and a value of 0 is returned.
   /// Similar the variable names are stored in the array strings.
   
   if(!jcalib){
      _DBG_<<"ERROR - GetConstants() called when jcalib not set!"<<endl;
      _DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }

   map <string, float> detparms;
   jcalib->Get(namepath, detparms);

   if((int)detparms.size()<*Nvals)
     *Nvals = (int)detparms.size();
   int i=0;
   for( map<string, float>::iterator ii=detparms.begin(); ii!=detparms.end(); ++ii){
     if (i<*Nvals){
       strcpy (strings[i].str, (*ii).first.c_str());
       vals[i++] = (*ii).second;
     }
   }
   return (int)detparms.size()>*Nvals; // return 0 if OK, 1 if not
}


//----------------
// GetColumn
//----------------
// Get a single column from the database by its key (string)
int GetColumn(const char* namepath, int *Nvals, float* vals, char *key_cstr){
 
   if(!jcalib){
      _DBG_<<"ERROR - GetColumn() called when jcalib not set!"<<endl;
      _DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }

   vector<map <string, float> >detparms;
   jcalib->Get(namepath, detparms);

   if (*Nvals!=int(detparms.size())){
     _DBG_ << "ERROR - Array size mismatch:  " << *Nvals << " != " 
	   << detparms.size() 
	   << " (ccdb)" << endl;
     _DBG_ << " for " << namepath <<endl;
     exit(-1); 
   }

   string key=key_cstr;

   unsigned int i=0;
   for (i=0;i<detparms.size();i++){     
     map<string, float> &row = detparms[i];
     vals[i]=row[key];
   }

   return 0;
}




//----------------
// GetArrayConstants
//----------------
int GetArrayConstants(const char* namepath, int *Nvals, float* vals, mystr_t *strings)
{
   /// C-callable routine for accessing calibration constants.
   /// The values specified by "namepath" will be read into the array
   /// "vals". The "vals" array should have enough memory allocated
   /// to hold *Nvals elements. If not, only the first *Nvals elements
   /// will be copied and a non-zero value returned. If the number
   /// of values in the database are less than *Nvals, then all values
   /// are copied, *Nvals is updated to reflect the number of valid
   /// elements in "vals", and a value of 0 is returned.
   /// Similar the variable names are stored in the array strings.
   
   if(!jcalib){
      _DBG_<<"ERROR - GetArrayConstants() called when jcalib not set!"<<endl;
      _DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
      _DBG_<<"ERROR - Exiting ..."<<endl;
      exit(-1);
   }

   vector<map <string, float> >detparms;
   jcalib->Get(namepath, detparms);

   unsigned int i=0;
   int j=0;
   for (i=0;i<detparms.size();i++){
     for( map<string, float>::iterator ii=detparms[i].begin(); ii!=detparms[i].end(); ++ii){
       if (j<*Nvals){
         strcpy (strings[j].str, (*ii).first.c_str());
         vals[j] = (*ii).second;
         j++;
       }
       else return 1;
     }
   }
   *Nvals=j;

   return 0; // return 0 if OK, 1 if not
}

//------------------
// GetMD5Geom
//------------------
const char* GetMD5Geom(void)
{
   // Get the MD5 checksum of the geometry that will be
   // used for the simulation. This will retrieve the
   // geometry checksum from either what has been statically
   // linked in, or dynamically, whichever is being used.

   // This is a little odd since the string originates
   // in a FORTRAN routine.
   static char md5[256];
   memset(md5, 0, 256);
   md5geom_wrapper_(md5);
   
   md5[32] = 0; // truncate string at 32 characters (FORTRAN adds a space)
   
   return md5;
}

extern "C" {

   #include <JANA/Calibrations/JCalibration.h>

   extern float coherent_peak_GeV[2];
   extern float endpoint_calib_GeV;
   extern float endpoint_energy_GeV;

   float get_endpoint_energy_(int &runno)
   {
      char dbname[] = "/PHOTON_BEAM/endpoint_energy";
      unsigned int ndata = 1;
      if (jcalib == 0)
         jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      if (GetCalib(dbname, &ndata, &endpoint_energy_GeV)) {
         fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
                 "failed to read photon beam endpoint energy",
                 "from ccdb, cannot continue.");
         exit (2);
      }
      return endpoint_energy_GeV;
   }
   
   float get_calib_energy_(int &runno)
   {
      char dbname[] = "/PHOTON_BEAM/hodoscope/endpoint_calib";
      unsigned int ndata = 1;
      if (jcalib == 0)
         jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      if (GetCalib(dbname, &ndata, &endpoint_calib_GeV)) {
         fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
                 "failed to read photon beam endpoint_calib",
                 "from ccdb, cannot continue.");
         exit (2);
      }
      return endpoint_calib_GeV;
   }
   
   float get_peak_energy_(int &runno)
   {
      char dbname[] = "/PHOTON_BEAM/coherent_energy";
      unsigned int ndata = 2;
      if (jcalib == 0)
         jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      if (GetCalib(dbname, &ndata, coherent_peak_GeV)) {
         fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
                 "failed to read coherent peak energy",
                 "from ccdb, cannot continue.");
         exit (2);
      }
      return coherent_peak_GeV[1];
   }
}
