// Author: David Lawrence Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.h
//
/// Processor for mcsmear
///

#ifndef _MYPROCESSOR_H_
#define _MYPROCESSOR_H_

#include <string>

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <fstream>
#include <HDDM/hddm_s.hpp>

#include "smear.h"
#include "mcsmear_config.h"

class MyProcessor:public JEventProcessor
{
   public:
   	  MyProcessor(mcsmear_config_t *in_config) {
   	  	 config = in_config;
   	  	 //smearer = NULL;
   	  }
   
      jerror_t init(void);                              ///< Called once at program start.
      jerror_t brun(JEventLoop *loop, int32_t runnumber);  ///< Called everytime a new run number is detected.
      jerror_t evnt(JEventLoop *loop, uint64_t eventnumber); ///< Called every event.
      jerror_t erun(void) {                             ///< Called everytime run number changes, provided brun has been called.
         return NOERROR;
      }
      jerror_t fini(void);                              ///< Called after last event of last event source has been processed.

      ofstream *ofs;
      hddm_s::ostream *fout; 
      unsigned long Nevents_written;

   private:
      int  HDDM_USE_COMPRESSION;
      bool HDDM_USE_INTEGRITY_CHECKS;
      
      mcsmear_config_t *config;
      //Smear *smearer;              //moved into thread_local storage, see MyProcessor.cc
};


#endif  // _MYPROCESSOR_H_
