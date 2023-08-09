C
C--- External beam energy distribution (alternative to the Bremsstrahlung beam) 
C--- The beam energy limits are defined elsewhere: EPH_LIM(2)
C
      INTEGER MXBMPN
      PARAMETER (MXBMPN=1000)
      COMMON/BG_BEAMEXT/ NBMEXT           ! number of points
     +                  ,IBMEXT           ! flag  =0 - equidistant points in EPH_LIM(1) - EPH_LIM(2) interval; 
C                                         !       =1 - the mean energies and the widths of the bin are set, will redefine EPH_LIM
     +                  ,EXT_LIM(2)       ! energy limits
     +                  ,RABMEXT(MXBMPN)  ! beam rate dN/dE
     +                  ,ENBMEXT(MXBMPN)  ! central energies of the bins (for IBMEXT>0) in GeV
C
      INTEGER            NBMEXT
     +                  ,IBMEXT
      REAL               EXT_LIM
     +                  ,RABMEXT
     +                  ,ENBMEXT
C