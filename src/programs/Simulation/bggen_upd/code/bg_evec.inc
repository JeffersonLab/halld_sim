C
C---  Simulated event
C
      INTEGER MXTRA
      PARAMETER (MXTRA=100)
      COMMON/EV_RECORD/ IEVENT,INDUM(2),IEVPROC
     +                 ,ITPIN(2,2),AMIN(2),PIN(3,2)
     +                 ,NTRA 
     +                 ,ITPTRA(6,MXTRA),AMTRA(MXTRA),PTRA(3,MXTRA)
      COMMON/EV_RECNAM/ CPROCNAM
C
      INTEGER IEVENT   ! event number
     +       ,INDUM    !  dummies (for later use)
     +       ,IEVPROC  ! the process number (=0 - PYTHIA)
     +       ,ITPIN    ! (1,i)=GEANT type, (2,i) - KF (LUND), k=1 - beam, =2 -target 
     +       ,NTRA     ! number of particles excluding the beam and the target
     +       ,ITPTRA   ! (1,k) track type (GEANT), (2-6,k) - LUND flags ; k=1,NTRA 
C                       (2,k) - KS code, 1-undecayed particle, 11 - decayed partice, others for various internal structures
C                       (3,k) KF (PDG) code
C                       (4,k) the entry number for the parent particle
C                       (5-6,k) the entry numbers for the first and last decay products
      REAL    AMIN     ! masses of the beam and the target
     +       ,PIN      ! (1-3,k) - 3-momenta, k=1 - beam, k=2 - target
     +       ,AMTRA    ! secondary particles' masses
     +       ,PTRA     ! 3-momenta
C
      CHARACTER CPROCNAM*48  ! the name of the process

 
