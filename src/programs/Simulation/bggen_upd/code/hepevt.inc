C
C HepEvt format
C used to interface with PHOTOS
C
C*******************************************************
C This is the "official" event common of the HEP community for MC events.
C  NEVHEP      = event number
C  NHEP        = number of particles ("lines")
C  ISTHEP(N)   = stability code; 0 if stable, 1 if decays.
C  IDHEP(I)    = particle ID code; see PDG '92, III.73
C  JMOHEP(1,I) = line number of mother particle
C  JMOHEP(2,I) = line number of last mother particle, if formed by string
C  JDAHEP(1,I) = line number of 1st daughter particle
C  JDAHEP(2,I) = line number of last daughter particle - MUST BE CONSECUTIVE!
C  PHEP(5,I)   = px,py,pz,E,m of particle
C  VHEP(4,I)   = x,y,z,L of production vertex
      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION  PHEP,VHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     *JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      SAVE /HEPEVT/
