/*
 * hitTag - registers hits for the tagger focal plane counters
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 * version 1.0 	-Richard Jones November 16, 2006
 * version 2.0 	-Richard Jones July 1, 2014
 * version 3.0  -A.S., October 16, 2019
 *
 * Programmer's Notes:
 * -------------------
 * 1) There is no tagger in the HDGeant simulation so no tagger hits are
 *    generated during tracking.  This hitTagger() function is called at
 *    event initialization time to register the tagged photon that is
 *    supposed to have caused the event. 
 * 2) Only microscope hits are produced in this version.
 * 3) In the simulation of physics events (external generator) with
 *    background enabled, pickTagger() produces a list of tagger hits
 *    that includes the original photon from the generator plus all
 *    of the background photons.  Note that this includes many photons
 *    that never reached the GlueX target because they were stopped
 *    at the collimator.
 *
 * update July 1, 2014 (version 2.0)
 * ---------------------------------
 * 1) Read the tagger channel energy bounds from the ccdb instead of
 *    hard-wiring them here.
 * 2) Add hits in both the fixed_array and microscope detectors.
 * 3) Fix the bug that forced the E value written into the hits structure
 *    to always contain the exact simulated beam photon energy, instead
 *    of the mean value for the hit tagger channel. Now only the mean
 *    photon energy for the hit channel is recorded.
 * 4) The recorded photon energy from the tagger is computed from the
 *    endpoint energy in the ccdb multiplied by the scaled_energy_range
 *    array values.
 *
 *
 *  update October 16, 2019 (version 3.0), A.S.
 * ---------------------------------
 * Modify calculation of the photon beam energy to account 
 * for the fact that the energy of bremsstrahlung electrons detected
 * by each tagger counter does not depend on the electron beam energy.
 * The photon beam energy E_gamma has to be computed as
 *
 *   E_gamma = R * E_endpoint_calib  +  DE,  where
 * 
 *   R is the fractional energy of the beam photon corresponding to the
 *  specific tagger counter. R is stored in the 'scaled_energy_range' 
 *  CCDB table. R is determined using calibration runs with the 
 *  specific electron beam energy, E_endpoint_calib
 *
 *   DE = Ebeam - E_endpoint_calib
 *
 *  Run-by-run dependent Ebeam is stored in the endpoint_energy CCDB table
 *  
 * 
 *  E_endpoint_calib energy is placed in the 
 *  /PHOTON_BEAM/hodoscope/endpoint_calib CCDB table
 *
 * for each specific run period. If the run (period) is not listed in the
 * table, the photon beam energy is computed as:
 *
 * 
 *  E_gamma = R * Ebeam 
 *
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <calibDB.h>

#define MICRO_TWO_HIT_RESOL     25.
#define MICRO_MAX_HITS          5000
#define FIXED_TWO_HIT_RESOL     25.
#define FIXED_MAX_HITS          5000
#define C_CM_PER_NS             29.9792458
#define TAG_T_MIN_NS            -200
#define TAG_T_MAX_NS            +200

float coherent_peak_GeV[2] = {0,0};
float endpoint_calib_GeV = 0;
float endpoint_energy_GeV = 0;
float micro_limits_Erange[2];
float hodo_limits_Erange[2];
static int micro_nchannels = 102;
float* micro_channel_Erange = 0;
static int hodo_nchannels = 274;
float* hodo_channel_Erange = 0;
binTree_t* microTree = 0;
binTree_t* hodoTree = 0;
static int microCount = 0;
static int hodoCount = 0;
static int printDone = 0;
static float beam_period = -1.0;

static int read_endpoint_calib = 0;


float get_reference_plane_();

/* register hits during event initialization (from gukine) */

void hitTagger (float xin[4], float xout[4],
                float pin[5], float pout[5], float dEsum,
                int track, int stack, int history)
{

   /* read beam_period from calibdb */
   if(beam_period < 0.0)
   {
          char dbname[] = "/PHOTON_BEAM/RF/beam_period::mc";
          unsigned int ndata = 1;
          if (GetCalib(dbname, &ndata, &beam_period)) {
        	 fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
        			 "failed to read RF period ",
        			 "from calibdb, cannot continue.");
        	 exit (2);
          }
   }

   int micro_chan;
   int hodo_chan;
   double Etag = 0;
   double E = pin[3];
   float ref_time_z_cm = get_reference_plane_();
   double t = xin[3]*1e9-(xin[2]-ref_time_z_cm)/C_CM_PER_NS;
   t = floor(t/beam_period+0.5)*beam_period;

   /* read tagger set endpoint energy from calibdb */
   if (endpoint_energy_GeV == 0) {
      char dbname[] = "/PHOTON_BEAM/endpoint_energy";
      unsigned int ndata = 1;
      if (GetCalib(dbname, &ndata, &endpoint_energy_GeV)) {
         fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
                 "failed to read photon beam endpoint energy",
                 "from calibdb, cannot continue.");
         exit (2);
      }
   }
 

   if (read_endpoint_calib == 0) {
      char dbname[] = "/PHOTON_BEAM/hodoscope/endpoint_calib";
      unsigned int ndata = 1;
      if (GetCalib(dbname, &ndata, &endpoint_calib_GeV)) {
        fprintf(stderr,"HDGeant error in hitTagger: failed to read endpoint_calib energy \n");
      }
      

      read_endpoint_calib = 1;

   }


   /* read microscope channel energy bounds from calibdb */
   if (micro_channel_Erange == 0) {

     int i;

      char dbname[] = "/PHOTON_BEAM/microscope/scaled_energy_range";
      /* table microscope/scaled_energy_range has 3 columns:
       *     column  xlow  xhigh
       * which are returned in an array like float[3][ncolumns]
       */
      int ndata = 3*micro_nchannels;
      mystr_t names[ndata];
      micro_channel_Erange = malloc(ndata*sizeof(float));
      if (GetArrayConstants(dbname, &ndata, micro_channel_Erange, names) ||
          ndata != 3*micro_nchannels)
      {
         fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
                 "failed to read microscope scaled_energy_range table",
                 "from calibdb, cannot continue.");
         exit (2);
      }
      else {

        micro_limits_Erange[0] = 0;
        micro_limits_Erange[1] = 1;
        for (i=0; i < micro_nchannels; ++i) {
          if (micro_limits_Erange[0] < micro_channel_Erange[3*i+1])
            micro_limits_Erange[0] = micro_channel_Erange[3*i+1];
          if (micro_limits_Erange[1] > micro_channel_Erange[3*i+2])
            micro_limits_Erange[1] = micro_channel_Erange[3*i+2];
        }
        
         
        if(endpoint_energy_GeV > 0 && endpoint_calib_GeV > 0) {

          printf(" Correct Beam Photon Energy  (TAGM) \n\n");
          
          double delta_E = endpoint_energy_GeV  -  endpoint_calib_GeV;
          
          for (i = 0; i < micro_nchannels; ++i) {
            micro_channel_Erange[3*i+1]  = micro_channel_Erange[3*i+1]*endpoint_calib_GeV  +  delta_E;
            micro_channel_Erange[3*i+2]  = micro_channel_Erange[3*i+2]*endpoint_calib_GeV  +  delta_E;
          }
          
          micro_limits_Erange[0]  =  micro_limits_Erange[0]*endpoint_calib_GeV  +  delta_E;
          micro_limits_Erange[1]  =  micro_limits_Erange[1]*endpoint_calib_GeV  +  delta_E;
          
        } else {
          
          for (i = 0; i < micro_nchannels; ++i) {
            
            micro_channel_Erange[3*i+1]  *=  endpoint_energy_GeV; 
            micro_channel_Erange[3*i+2]  *=  endpoint_energy_GeV;	  
          }
          
          micro_limits_Erange[0] *= endpoint_energy_GeV;
          micro_limits_Erange[1] *= endpoint_energy_GeV;	
        }
               	
      }
   }
 

   /* read hodoscope channel energy bounds from calibdb */
   if (hodo_channel_Erange == 0) {

     int i;

      char dbname[] = "/PHOTON_BEAM/hodoscope/scaled_energy_range";
      /* table hodoscope/scaled_energy_range has 3 columns:
       *     counter  xlow  xhigh
       * which are returned in an array like float[3][ncolumns]
       */
      int ndata = 3*hodo_nchannels;
      mystr_t names[ndata];

      hodo_channel_Erange = malloc(ndata*sizeof(float));

      if (GetArrayConstants(dbname, &ndata, hodo_channel_Erange, names) ||
          ndata != 3*hodo_nchannels)
        {
          fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
        	  "failed to read hodoscope scaled_energy_range table",
        	  "from calibdb, cannot continue.");
          exit (2);
        }
      else {
        
        hodo_limits_Erange[0] = 0;
        hodo_limits_Erange[1] = 1;
        for (i=0; i < hodo_nchannels; ++i) {
          if (hodo_limits_Erange[0] < hodo_channel_Erange[3*i+1])
            hodo_limits_Erange[0] = hodo_channel_Erange[3*i+1];
          if (hodo_limits_Erange[1] > hodo_channel_Erange[3*i+2])
            hodo_limits_Erange[1] = hodo_channel_Erange[3*i+2];
        }
      }


      if(endpoint_energy_GeV > 0 && endpoint_calib_GeV > 0) {

        printf(" Correct Beam Photon Energy  (TAGH) \n\n");
        
        double delta_E = endpoint_energy_GeV  -  endpoint_calib_GeV;
        
        for (i = 0; i < hodo_nchannels; ++i) {
          hodo_channel_Erange[3*i + 1] = hodo_channel_Erange[3*i+1]*endpoint_calib_GeV  +  delta_E;
          hodo_channel_Erange[3*i + 2] = hodo_channel_Erange[3*i+2]*endpoint_calib_GeV  +  delta_E;
        }
        
        hodo_limits_Erange[0] = hodo_limits_Erange[0]*endpoint_calib_GeV + delta_E;
        hodo_limits_Erange[1] = hodo_limits_Erange[1]*endpoint_calib_GeV + delta_E;
        
      } else {	
        
        for (i = 0; i < hodo_nchannels; ++i) {
          hodo_channel_Erange[3*i+1] *= endpoint_energy_GeV;
          hodo_channel_Erange[3*i+2] *= endpoint_energy_GeV;
        }
        
        hodo_limits_Erange[0] *= endpoint_energy_GeV;
        hodo_limits_Erange[1] *= endpoint_energy_GeV;
        
      }          
      
   }
   



   if (printDone == 0) {
      fprintf(stderr,"TAGGER: ALL parameters loaded from Data Base\n");
      printDone = 1;
   }

   /* look up hit tagger channel, if any */
   hodo_chan = -1;
   micro_chan = -1;
   if (E < micro_limits_Erange[0] && E > micro_limits_Erange[1]) {
      int i;
      for (i=0; i < micro_nchannels; ++i) {
         if ( E < micro_channel_Erange[3*i+1] &&
              E > micro_channel_Erange[3*i+2] )
         {
            Etag = (micro_channel_Erange[3*i+1] + 
                    micro_channel_Erange[3*i+2]) / 2;
            micro_chan = micro_channel_Erange[3*i];
            break;
         }
      }
   }
   else if (E < hodo_limits_Erange[0] && E > hodo_limits_Erange[1]) {
      int i;
      for (i=0; i < hodo_nchannels; ++i) {
         if ( E < hodo_channel_Erange[3*i+1] &&
              E > hodo_channel_Erange[3*i+2] )
         {
            Etag = (hodo_channel_Erange[3*i+1] +
                    hodo_channel_Erange[3*i+2]) / 2;
            hodo_chan = hodo_channel_Erange[3*i];
            break;
         }
      }
   }

   /* post the hit to the microscope hits tree, mark channel as hit */

   if (micro_chan > -1) {
      int nhit;
      s_TaggerTruthHits_t* hits;
      int mark = micro_chan + 1000;
      void** twig = getTwig(&microTree, mark);
      if (*twig == 0)
      {
         s_Tagger_t* tag = *twig = make_s_Tagger();
         s_MicroChannels_t* channels = make_s_MicroChannels(1);
         hits = make_s_TaggerTruthHits(MICRO_MAX_HITS);
         hits->mult = 0;
         channels->in[0].taggerTruthHits = hits;
         channels->in[0].column = micro_chan;
         channels->in[0].row = 0;
         channels->in[0].E = Etag;
         channels->mult = 1;
         tag->microChannels = channels;
         microCount++;
      }
      else
      {
         s_Tagger_t* tag = *twig;
         hits = tag->microChannels->in[0].taggerTruthHits;
      }
   
      if (hits != HDDM_NULL)
      {
         for (nhit = 0; nhit < hits->mult; nhit++)
         {
            if (fabs(hits->in[nhit].t - t) < MICRO_TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < hits->mult)         /* ignore second hit */
         {
         }
         else if (nhit < MICRO_MAX_HITS)   /* create new hit */
         {
            hits->in[nhit].bg = track;
            hits->in[nhit].t = t;
            hits->in[nhit].E = E;
            hits->in[nhit].dE += 3.5e-3; // GeV in SciFi
            hits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitTagger: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",
                    MICRO_MAX_HITS);
         }
      }
   }

   /* post the hit to the hodoscope hits tree, mark channel as hit */

   if (hodo_chan > -1) {
      int nhit;
      s_TaggerTruthHits_t* hits;
      int mark = hodo_chan + 1000;
      void** twig = getTwig(&hodoTree, mark);
      if (*twig == 0)
      {
         s_Tagger_t* tag = *twig = make_s_Tagger();
         s_HodoChannels_t* channels = make_s_HodoChannels(1);
         hits = make_s_TaggerTruthHits(FIXED_MAX_HITS);
         hits->mult = 0;
         channels->in[0].taggerTruthHits = hits;
         channels->in[0].counterId = hodo_chan;
         channels->in[0].E = Etag;
         channels->mult = 1;
         tag->hodoChannels = channels;
         hodoCount++;
      }
      else
      {
         s_Tagger_t* tag = *twig;
         hits = tag->hodoChannels->in[0].taggerTruthHits;
      }
   
      if (hits != HDDM_NULL)
      {
         for (nhit = 0; nhit < hits->mult; nhit++)
         {
            if (fabs(hits->in[nhit].t - t) < FIXED_TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < hits->mult)         /* ignore second hit */
         {
         }
         else if (nhit < FIXED_MAX_HITS)   /* create new hit */
         {
            hits->in[nhit].bg = track;
            hits->in[nhit].t = t;
            hits->in[nhit].E = E;
            hits->in[nhit].dE += 5.5e-4; // GeV in hodo scint.
            hits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitTagger: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",
                    FIXED_MAX_HITS);
         }
      }
   }
}

/* entry point from fortran */

void hittagger_ (float* xin, float* xout,
                 float* pin, float* pout, float* dEsum,
                 int* track, int* stack, int* history)
{
   hitTagger(xin,xout,pin,pout,*dEsum,*track,*stack,*history);
}


/* pick and package the hits for shipping */

s_Tagger_t* pickTagger ()
{
   s_Tagger_t* box;
   s_Tagger_t* item;

   if (microCount == 0 && hodoCount == 0)
   {
      return HDDM_NULL;
   }

   box = make_s_Tagger();

   box->microChannels = make_s_MicroChannels(microCount);
   while ((item = (s_Tagger_t*) pickTwig(&microTree)))
   {
      s_MicroChannels_t* channels = item->microChannels;
      int channel;
      for (channel=0; channel < channels->mult; ++channel)
      {
         s_TaggerTruthHits_t* hits = channels->in[channel].taggerTruthHits;

         /* constraint t values to lie within time range */
         int i;
         int iok=0;
         for (iok=i=0; i < hits->mult; i++)
         {
            if ((hits->in[i].t >= TAG_T_MIN_NS) &&
                (hits->in[i].t <= TAG_T_MAX_NS))
            {
               if (iok < i)
               {
                  hits->in[iok] = hits->in[i];
               }
               ++iok;
            }
         }
         if (iok)
         {
            hits->mult = iok;
            int m = box->microChannels->mult++;
            box->microChannels->in[m] = channels->in[0];
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (channels != HDDM_NULL)
      {
         FREE(channels);
      }
      FREE(item);
   }

   box->hodoChannels = make_s_HodoChannels(hodoCount);
   while ((item = (s_Tagger_t*) pickTwig(&hodoTree)))
   {
      s_HodoChannels_t* channels = item->hodoChannels;
      int channel;
      for (channel=0; channel < channels->mult; ++channel)
      {
         s_TaggerTruthHits_t* hits = channels->in[channel].taggerTruthHits;

         /* constraint t values to lie within time range */
         int i;
         int iok=0;
         for (iok=i=0; i < hits->mult; i++)
         {
            if ((hits->in[i].t >= TAG_T_MIN_NS) &&
                (hits->in[i].t <= TAG_T_MAX_NS))
            {
               if (iok < i)
               {
                  hits->in[iok] = hits->in[i];
               }
               ++iok;
            }
         }
         if (iok)
         {
            hits->mult = iok;
            int m = box->hodoChannels->mult++;
            box->hodoChannels->in[m] = channels->in[0];
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (channels != HDDM_NULL)
      {
         FREE(channels);
      }
      FREE(item);
   }

   microCount = 0;
   hodoCount = 0;

   if ((box->microChannels != HDDM_NULL) &&
       (box->microChannels->mult == 0))
   {
      FREE(box->microChannels);
      box->microChannels = HDDM_NULL;
   }
   if ((box->hodoChannels != HDDM_NULL) &&
       (box->hodoChannels->mult == 0))
   {
      FREE(box->hodoChannels);
      box->hodoChannels = HDDM_NULL;
   }
   if (box->microChannels->mult == 0 &&
       box->hodoChannels->mult == 0)
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
