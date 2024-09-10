/*
 * hitCCal - registers hits for Compton calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 *      version 1.1     -A. Somov  
 *      Add interface with CCDB,  change module length and attenuation length
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <gid_map.h>

#include "calibDB.h"

extern s_HDDM_t* thisInputEvent;


static float   ATTEN_LENGTH     =   200.;     //effective attenuation length in PbWO
static float   C_EFFECTIVE      =   13.;      //effective speed of light in PbWO cm/ns
static float   WIDTH_OF_BLOCK   =   2.;       //cm
static float   LENGTH_OF_BLOCK  =   20.;      //cm
static float   TWO_HIT_RESOL    =   75.;      //ns
static int     MAX_HITS         =   100;
static float   THRESH_MEV       =   5.;



binTree_t* ComptonCalTree = 0;
static int blockCount = 0;
static int showerCount = 0;

static int initialized = 0;

/* register hits during tracking (from gustep) */

void hitComptonEMcal (float xin[4], float xout[4],
                    float pin[5], float pout[5], float dEsum,
                    int track, int stack, int history, int ipart)
{
   float x[3], t;
   float xccal[3];


  if (!initialized){
     
     mystr_t strings[50];
     float  values[50];
     int  nvalues = 50;
     int status = GetConstants("CCAL/ccal_parms", &nvalues, values, strings);

     if (!status) {
       
       int ncounter = 0;
       int i;
       for ( i = 0; i < (int)nvalues; i++){

	 if (!strcmp(strings[i].str,"CCAL_ATTEN_LENGTH")) {
	   ATTEN_LENGTH  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_C_EFFECTIVE")) {
	   C_EFFECTIVE  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_WIDTH_OF_BLOCK")) {
	   WIDTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_LENGTH_OF_BLOCK")) {
	   LENGTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_TWO_HIT_RESOL")) {
	   TWO_HIT_RESOL  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_MAX_HITS")) {
	   MAX_HITS  = (int)values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"CCAL_THRESH_MEV")) {
	   THRESH_MEV  = values[i];
	   ncounter++;
	 }

       }


       const int nparams = 7;

       if (ncounter == nparams){
	 printf("CCAL: ALL parameters loaded from Data Base\n");
       } else if (ncounter < nparams){
         printf("CCAL: NOT ALL necessary parameters found in Data Base %d out of %d\n",ncounter,nparams);
       } else {
	 printf("CCAL: SOME parameters found more than once in Data Base\n");
       }
     } else {
       printf("CCAL: Cannot Load Parameters from Data Base. Use default values. \n");
     }

     initialized = 1;

   }


   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xccal,"CCAL");

   /* post the hit to the truth tree */

   int itrack = (stack == 0)? gidGetId(track) : -1;

   if ((history == 0) && (pin[3] > THRESH_MEV/1e3))
   {
      s_CcalTruthShowers_t* showers;
      int mark = (1<<30) + showerCount;
      void** twig = getTwig(&ComptonCalTree, mark);
      if (*twig == 0)
      {
         s_ComptonEMcal_t* cal = *twig = make_s_ComptonEMcal();
         cal->ccalTruthShowers = showers = make_s_CcalTruthShowers(1);
         int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         showers->in[0].primary = (track <= a && stack == 0);
         showers->in[0].track = track;
         showers->in[0].t = xin[3]*1e9;
         showers->in[0].x = xin[0];
         showers->in[0].y = xin[1];
         showers->in[0].z = xin[2];
         showers->in[0].px = pin[0]*pin[4];
         showers->in[0].py = pin[1]*pin[4];
         showers->in[0].pz = pin[2]*pin[4];
         showers->in[0].E = pin[3];
         showers->in[0].ptype = ipart;
         showers->in[0].trackID = make_s_TrackID();
         showers->in[0].trackID->itrack = itrack;
         showers->mult = 1;
         showerCount++;
      }
   }

   /* post the hit to the hits tree, mark block as hit */

   if (dEsum > 0)
   {
      int nhit;
      s_CcalTruthHits_t* hits;
      int row = getrow_wrapper_();
      int column = getcolumn_wrapper_();
      
      float dist = 0.5*LENGTH_OF_BLOCK-xccal[2];
      float dEcorr = dEsum * exp(-dist/ATTEN_LENGTH);
      float tcorr = t + dist/C_EFFECTIVE;
      int mark = ((row+1)<<16) + (column+1);
      void** twig = getTwig(&ComptonCalTree, mark);
      if (*twig == 0)
      {
         s_ComptonEMcal_t* cal = *twig = make_s_ComptonEMcal();
         s_CcalBlocks_t* blocks = make_s_CcalBlocks(1);
         blocks->mult = 1;
         blocks->in[0].row = row;
         blocks->in[0].column = column;
         blocks->in[0].ccalTruthHits = hits = make_s_CcalTruthHits(MAX_HITS);
         cal->ccalBlocks = blocks;
         blockCount++;
      }
      else
      {
         s_ComptonEMcal_t* cal = *twig;
         hits = cal->ccalBlocks->in[0].ccalTruthHits;
      }

      for (nhit = 0; nhit < hits->mult; nhit++)
      {
         if (fabs(hits->in[nhit].t - tcorr) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < hits->mult)		/* merge with former hit */
      {
			/* unclear if the intent here was to add dEcorr to hits->in[nhit].E */
			/* in the numerator as well as denominator. This caused a compiler  */
			/* warning so I chose for it not to. (I'm pretty sure that's right) */
			/*                                    10/28/2015 DL                 */
			/*
         hits->in[nhit].t =
                       (hits->in[nhit].t * hits->in[nhit].E + tcorr*dEcorr)
                     / (hits->in[nhit].E += dEcorr);
		   */
			hits->in[nhit].t =
                       (hits->in[nhit].t * hits->in[nhit].E + tcorr*dEcorr)
                     / (hits->in[nhit].E + dEcorr);
			hits->in[nhit].E += dEcorr;
      }
      else if (nhit < MAX_HITS)         /* create new hit */
      {
         hits->in[nhit].t = tcorr;
         hits->in[nhit].E = dEcorr;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitComptonEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         //exit(2);
      }
   }
}

/* entry point from fortran */

void hitcomptonemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum,
                      int* track, int* stack, int* history, int* ipart)
{
   hitComptonEMcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


/* pick and package the hits for shipping */

s_ComptonEMcal_t* pickComptonEMcal ()
{
   s_ComptonEMcal_t* box;
   s_ComptonEMcal_t* item;

   if ((blockCount == 0) && (showerCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_ComptonEMcal();
   box->ccalBlocks = make_s_CcalBlocks(blockCount);
   box->ccalTruthShowers = make_s_CcalTruthShowers(showerCount);
   while ((item = (s_ComptonEMcal_t*) pickTwig(&ComptonCalTree)))
   {
      s_CcalBlocks_t* blocks = item->ccalBlocks;
      int block;
      s_CcalTruthShowers_t* showers = item->ccalTruthShowers;
      int shower;
      for (block=0; block < blocks->mult; ++block)
      {
         s_CcalTruthHits_t* hits = blocks->in[block].ccalTruthHits;

         if (hits)
         {
            int m = box->ccalBlocks->mult;

         /* compress out the hits below threshold */
            int i,iok;
            for (iok=i=0; i < hits->mult; i++)
            {
               if (hits->in[i].E > THRESH_MEV/1e3)
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
               box->ccalBlocks->in[m] = blocks->in[block];
               box->ccalBlocks->mult++;
            }
            else if (hits != HDDM_NULL)
            {
               FREE(hits);
            }
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }

      for (shower=0; shower < showers->mult; ++shower)
      {
         int m = box->ccalTruthShowers->mult++;
         box->ccalTruthShowers->in[m] = showers->in[shower];
      }
      if (blocks != HDDM_NULL)
      {
         FREE(blocks);
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }
      FREE(item);
   }

   blockCount = showerCount = 0;

   if ((box->ccalBlocks != HDDM_NULL) &&
       (box->ccalBlocks->mult == 0))
   {
      FREE(box->ccalBlocks);
      box->ccalBlocks = HDDM_NULL;
   }
   if ((box->ccalTruthShowers != HDDM_NULL) &&
       (box->ccalTruthShowers->mult == 0))
   {
      FREE(box->ccalTruthShowers);
      box->ccalTruthShowers = HDDM_NULL;
   }
   if ((box->ccalBlocks->mult == 0) &&
       (box->ccalTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
