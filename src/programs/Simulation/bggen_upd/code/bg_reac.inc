C
C---  A binary reaction definition: gamma+A-->B+C(recoil) or 1+2-->3+4
C---     t=(p3-p1)^2
C
      INTEGER MXPNTR,MXTRNG
      PARAMETER (MXPNTR=3000,MXTRNG=5)
      COMMON/BG_REAC/  IPREAC(2)
     +                ,TLREAC(2)
     +                ,TPAREAC(11,MXTRNG)
     +                ,ITDRNG(MXTRNG),ITDSING,TSLREAC
     +                ,ELREAC(2),NPXREAC
     +                ,XSREAC(MXPNTR),ESREAC(MXPNTR)
C
      INTEGER IPREAC      ! (1:2) - GEANT types of the secondary particles (=0 - photon) (F)
     +       ,NPXREAC     ! number of points in the array of energies/cross_sections (F)
     +       ,ITDRNG      ! (irng)>0 if t-dist range irng is active (calculated), = 2 if a simple single exponential
     +       ,ITDSING     ! = irng of the "single" exponential, if exists, or 0
C
      REAL    ELREAC      ! (1:2) min,max of the beam energy (F)
     +       ,ESREAC      ! beam energies ELREAC(1)+(ELREAC(2)-ELREAC(1))/(NPXREAC-1)*(i-1) 
     +       ,XSREAC      ! (1:NPXREAC) cross sections in the points ESREAC (F)
     +       ,TLREAC      ! (1:2) min max of the t range (t is negative!) (F)
     +       ,TPAREAC     ! (1,1) = EMX1 for beam energies <EMX1 range for the irt-th t-distribution (F)
C                         ! (2:11,1)  parameters for the t distribution function a1,a2,a3,a4,...a10 (F)
C                         ! (1,2) = EMX2  (EMX1<EMX2 for this range to be used) .... parameters ....
C                         ! ..........
C                         ! (1,5) = EMX5 .....parameters .... 
     +       ,TSLREAC     ! the slope of the single exponential if ITDSIGN>0
C
C      t-slope 1/(GeV/C)**2, ds/dt~a1*(exp(a2*t)+exp(a3*t)*a4)/(1-a5*t)**4
C
