/*
 * hitECal - registers hits for Compton calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *      Initial implementation version 0.1   A.S.  12/18/2023   
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



binTree_t* CrystalCalTree = 0;
static int blockCount = 0;
static int showerCount = 0;

static int initialized = 0;

/* register hits during tracking (from gustep) */

void hitCrystalEcal (float xin[4], float xout[4],
                    float pin[5], float pout[5], float dEsum,
                    int track, int stack, int history, int ipart)
{
   float x[3], t;
   float xecal[3];


  if (!initialized){
     
     mystr_t strings[50];
     float  values[50];
     int  nvalues = 50;


#if 1

     int status = GetConstants("ECAL/ecal_parms", &nvalues, values, strings);

     //     printf(" DB status = %d \n",status);

     if (!status) {
       
       int ncounter = 0;
       int i;
       for ( i = 0; i < (int)nvalues; i++){

	 if (!strcmp(strings[i].str,"ECAL_ATTEN_LENGTH")) {
	   ATTEN_LENGTH  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_C_EFFECTIVE")) {
	   C_EFFECTIVE  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_WIDTH_OF_BLOCK")) {
	   WIDTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_LENGTH_OF_BLOCK")) {
	   LENGTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_TWO_HIT_RESOL")) {
	   TWO_HIT_RESOL  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_MAX_HITS")) {
	   MAX_HITS  = (int)values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"ECAL_THRESH_MEV")) {
	   THRESH_MEV  = values[i];
	   ncounter++;
	 }

       }

       printf("ATTEN_LENGTH =  %f  C_EFFECTIVE = %f  WIDTH_OF_BLOCK  = %f  LENGTH_OF_BLOCK = %f  ECAL_TWO_HIT_RESOL = %f  ECAL_MAX_HITS = %d  ECAL_THRESH_MEV = %f  \n", ATTEN_LENGTH, C_EFFECTIVE, WIDTH_OF_BLOCK, LENGTH_OF_BLOCK, TWO_HIT_RESOL, MAX_HITS, THRESH_MEV);

       const int nparams = 7;

       if (ncounter == nparams){
	 printf("ECAL: ALL parameters loaded from Data Base\n");
       } else if (ncounter < nparams){
         printf("ECAL: NOT ALL necessary parameters found in Data Base %d out of %d\n",ncounter,nparams);
       } else {
	 printf("ECAL: SOME parameters found more than once in Data Base\n");
       }
     } else {
       printf("ECAL: Cannot Load Parameters from Data Base. Use default values. \n");
     }

#endif

     initialized = 1;

   }


   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xecal,"ECAL");

   /* post the hit to the truth tree */

   int itrack = (stack == 0)? gidGetId(track) : -1;

   if ((history == 0) && (pin[3] > THRESH_MEV/1e3))
   {
      s_EcalTruthShowers_t* showers;
      int mark = (1<<30) + showerCount;
      void** twig = getTwig(&CrystalCalTree, mark);
      if (*twig == 0)
      {
         s_CrystalEcal_t* cal = *twig = make_s_CrystalEcal();
         cal->ecalTruthShowers = showers = make_s_EcalTruthShowers(1);
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
      s_EcalTruthHits_t* hits;
      int row = getrow_wrapper_();
      int column = getcolumn_wrapper_();
      
      float dist = 0.5*LENGTH_OF_BLOCK-xecal[2];
      float dEcorr = dEsum * exp(-dist/ATTEN_LENGTH);
      float tcorr = t + dist/C_EFFECTIVE;
      int mark = ((row+1)<<16) + (column+1);
      void** twig = getTwig(&CrystalCalTree, mark);
      if (*twig == 0)
      {
         s_CrystalEcal_t* cal   = *twig = make_s_CrystalEcal();
         s_EcalBlocks_t* blocks = make_s_EcalBlocks(1);
         blocks->mult = 1;
         blocks->in[0].row = row;
         blocks->in[0].column = column;
         blocks->in[0].ecalTruthHits = hits = make_s_EcalTruthHits(MAX_HITS);
         cal->ecalBlocks = blocks;
         blockCount++;
      }
      else
      {
         s_CrystalEcal_t*  cal = *twig;
         hits = cal->ecalBlocks->in[0].ecalTruthHits;
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
         fprintf(stderr,"HDGeant error in hitCrystalEcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitcrystalecal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum,
                      int* track, int* stack, int* history, int* ipart)
{
   hitCrystalEcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


/* pick and package the hits for shipping */

s_CrystalEcal_t* pickCrystalEcal ()
{
   s_CrystalEcal_t* box;
   s_CrystalEcal_t* item;

   if ((blockCount == 0) && (showerCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_CrystalEcal();
   box->ecalBlocks = make_s_EcalBlocks(blockCount);
   box->ecalTruthShowers = make_s_EcalTruthShowers(showerCount);
   while ((item = (s_CrystalEcal_t*) pickTwig(&CrystalCalTree)))
   {
      s_EcalBlocks_t* blocks = item->ecalBlocks;
      int block;
      s_EcalTruthShowers_t* showers = item->ecalTruthShowers;
      int shower;
      for (block=0; block < blocks->mult; ++block)
      {
         s_EcalTruthHits_t* hits = blocks->in[block].ecalTruthHits;

         if (hits)
         {
            int m = box->ecalBlocks->mult;

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
               box->ecalBlocks->in[m] = blocks->in[block];
               box->ecalBlocks->mult++;
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
         int m = box->ecalTruthShowers->mult++;
         box->ecalTruthShowers->in[m] = showers->in[shower];
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

   if ((box->ecalBlocks != HDDM_NULL) &&
       (box->ecalBlocks->mult == 0))
   {
      FREE(box->ecalBlocks);
      box->ecalBlocks = HDDM_NULL;
   }
   if ((box->ecalTruthShowers != HDDM_NULL) &&
       (box->ecalTruthShowers->mult == 0))
   {
      FREE(box->ecalTruthShowers);
      box->ecalTruthShowers = HDDM_NULL;
   }
   if ((box->ecalBlocks->mult == 0) &&
       (box->ecalTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
