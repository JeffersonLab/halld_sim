C
C---  Run statistics
C
      INTEGER MXSTAT
      PARAMETER (MXSTAT=10)
      COMMON/BG_STATC/ NEVSTAT(MXSTAT)
C
      INTEGER NEVSTAT  ! (1) - number of initial events simulated
     +                 ! (2) - successfully simulated events
     +                 ! (3) - selected events  
     +                 ! (4) - written HDDM events  
     +                 ! (5) - written NTUPLE events  
     +                 ! (6) - total number of beam particles simulated 
