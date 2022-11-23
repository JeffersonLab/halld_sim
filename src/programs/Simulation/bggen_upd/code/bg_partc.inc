C
C---  Particle masses, widths, and the decay properties for locally processed decays (ISIMUL>0)  
C---  A link to PYTHIA particle codes
C
      INTEGER MXPART,MXPDEF,MXPKC
      PARAMETER (MXPART=500,MXPDEF=200,MXPKC=1000)
      COMMON/BG_PARTC/ NPADEF,IDGDEF(MXPDEF)
     +                ,KF_PART(MXPART)
     +                ,KC_PART(MXPART)
     +                ,LD_PART(MXPART)
     +                ,AM_PART(MXPART),WG_PART(MXPART)
     +                ,KD_PART(4,MXPART)
     +                ,ND_PART(MXPART)
     +                ,IDG_KC(-MXPKC:MXPKC)
     +                ,KF_KC(-MXPKC:MXPKC)
      COMMON/BG_PARTCC/ CH_PART(MXPART)
C     +                 ,CH_KC(-MXPKC:MXPKC)
C
      INTEGER  NPADEF   ! nunber of defined particles, ip=1,NPADEF
     +        ,IDGDEF   ! (ip) = idgea  particle ID (GEANT id)
      INTEGER  KF_PART  ! (idgea) = KF (PYTHIA and PDG particle code) should .ne.0
     +        ,KC_PART  ! (idgea) = KC PYTHIA internal code (if PYTHIA is initialized - ISIMUL=0),
     +        ,LD_PART  ! (idgea) = 0 no decay of this KF in PYTHIA (for particles transferred to GEANT)
     +        ,KD_PART  ! (1-3,idgea) idgea_p of the decay products (one decay/idgea allowed) for the local decay simulation
C                       ! (4,idgea) =0 -  decay uniform in theta, =1 - like rho (sin**2), =2 - J/Psi-type
     +        ,ND_PART  ! (idgea)      - number of products (calculated from KD_PART)
     +        ,IDG_KC   ! (kc)  = ID Geant for the PYTHIA code kc, if PYTHIA is initialized
     +        ,KF_KC    ! (kc)  = KF PYTHIA codes
C
      REAL     AM_PART  !  (idgea) - particle mass GeV/c^2
     +        ,WG_PART  !  (idgea) - full widths
      CHARACTER CH_PART*16  ! (idgea) particle name
C      CHARACTER CH_KC*16    ! (kc) particle name from PYTHIA
C
      COMMON/BG_STAT/  MULT_KC(-MXPKC:MXPKC)
     +                ,NODC_KC(-MXPKC:MXPKC)
C
      INTEGER MULT_KC  ! (kc) - full multiplicity per run for the particle KC in PYTHIA
     +       ,NODC_KC  ! (kc) - number of not decayed (fragmented) particles among MULT_KC
C


