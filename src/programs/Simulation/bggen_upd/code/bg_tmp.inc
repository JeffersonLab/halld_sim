C
      COMMON/BG_PARTMP/ K_PTMP(128)             ! temporary storage for an FFREAD command PARTICLE
      INTEGER           K_PTMP
      REAL              R_PTMP(128)
      EQUIVALENCE      (K_PTMP(1),R_PTMP(1))
C
C--------Temporary storage of PYTHIA modification command
C
      INTEGER MXPYMOD               ! max number of PYTHIA modifications (command)
      PARAMETER (MXPYMOD=500)
      COMMON/BG_PYTMPN/ NPYMOD          ! number of actual modifications (FF commands)
     +                 ,KPYTMP(8)       ! temporary storage for an FF key *PYTMOD                
      INTEGER           NPYMOD
     +                 ,KPYTMP
      COMMON/BG_PYTMP/  CPYMOD(MXPYMOD) ! array of commands
      CHARACTER  CPYMOD*32
