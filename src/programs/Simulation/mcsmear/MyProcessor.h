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
   
      void Init() override;
      void BeginRun(const std::shared_ptr<const JEvent>& event) override;
      void Process(const std::shared_ptr<const JEvent>& event) override;
      void EndRun() override  {
         return; //NOERROR;
      }
      void Finish() override;
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
