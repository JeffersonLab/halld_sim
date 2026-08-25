C---  COMMON for the spectral function data
C---  2023/11/16
C
      INTEGER MXSFB
      PARAMETER (MXSFB=1000)
      COMMON/SF_DATA/   SF_A                 ! A mass of the nucleus
     +                 ,NE_SF                ! number of bins in E
     +                 ,NP_SF                ! number of bins in P
     +                 ,SF_PROB(MXSFB,MXSFB) ! (iE,iP) probability of the bin: = SF*p**2*4pi*dP*dE
     +                 ,SF_ECEN(MXSFB)       ! (i) - E center of the bin i
     +                 ,SF_PCEN(MXSFB)       ! (i) - P center of the bin i
     +                 ,SF_EWID(MXSFB)       ! (i) - E width  of the bin i
     +                 ,SF_PWID(MXSFB)       ! (i) - P width  of the bin i
C
      INTEGER           NE_SF
     +                 ,NP_SF
      REAL              SF_A
     +                 ,SF_PROB
     +                 ,SF_ECEN
     +                 ,SF_PCEN
     +                 ,SF_EWID
     +                 ,SF_PWID


	
