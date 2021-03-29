*
* $Id: gutrak.F,v 1.1.1.1 1995/10/24 10:21:46 cernlib Exp $
*
* $Log: gutrak.F,v $
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUTRAK
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to control tracking of one track            *
C.    *                                                                *
C.    *    ==>Called by : GTREVE                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctrak.inc"
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      integer saved_istak
      common /gutrak_saved_istak/saved_istak
      if (ISTAK.eq.0) then
         saved_istak = IQ(JSTAK+1) + 1
      endif
*
      CALL GTRACK
*                                                             END GUTRAK
      END
*
* Returns the stack index of a primary track
*
      function istak_of_primary()
      integer istak_of_primary
      integer saved_istak
      common /gutrak_saved_istak/saved_istak
      istak_of_primary = saved_istak
      end
