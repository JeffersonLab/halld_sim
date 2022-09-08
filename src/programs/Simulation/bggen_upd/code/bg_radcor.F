      SUBROUTINE BG_RADCOR
C
C---  Radiative correction to particle decays using PHOTOS
C
      IMPLICIT NONE
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_evec.inc'
      INCLUDE 'hepevt.inc'
C
      INTEGER i,j,idgea,ir,ih,ifind
      INTEGER nrc,igrc(MXTRA),kfrc(MXTRA)  ! for PHOTOS - store the relation between Id_geant and KF
      REAL qq
C
C     ------------------------------------------------------------------
C
      IF(NTRA.LE.0.OR.NRADCOR.LE.0) GO TO 999 
C
C---      Any particle in the record to apply rad correction?
C
      ifind=0
      DO ir=1,NRADCOR
         idgea=IRADCORR(ir)
         DO i=1,NTRA
            IF(ITPTRA(1,i).EQ.idgea) THEN
               ifind=1
               GO TO 10
            ENDIF
         ENDDO
      ENDDO
 10   IF(ifind.EQ.0) GO TO 999
C
C---        Store the IDG - KF relation for this record
C
      nrc=NTRA
      DO j=1,nrc
         igrc(j)=ITPTRA(1,j)
         kfrc(j)=ITPTRA(3,j)
      ENDDO
      nrc=nrc+1                 ! add a photon
      igrc(nrc)=1
      kfrc(nrc)=22
C
C---    Fill HepEvt block as input to PHOTOS
C
      NEVHEP = IEVENT
      NHEP = NTRA
      DO i=1,NHEP
         ISTHEP(i) = ITPTRA(2,i)
C         IF(ISTHEP(i).EQ.10) THEN
C            ISTHEP(i) = 2
C         ENDIF
         IDHEP(i) = ITPTRA(3,i)
         JMOHEP(1,i) = ITPTRA(4,i)
         JMOHEP(2,i) = 0
         JDAHEP(1,i) = ITPTRA(5,i)
         JDAHEP(2,i) = ITPTRA(6,i)
         qq=0.
         DO j=1,3
            PHEP(j,i) = PTRA(j,i)
            qq=qq+PTRA(j,i)**2
         ENDDO
         PHEP(4,i) = SQRT(AMTRA(i)**2+qq)
         PHEP(5,i) = AMTRA(i)
         DO j=1,4
            VHEP(j,i) = 0.
         ENDDO
      ENDDO
C      IF(IEVENT.LT.NPRIEV) THEN
C         WRITE(6,*) ' ==== Before PHOTOS '
C         CALL BG_PRI(1)
C      ENDIF
C
C---  Add radiation corrections
C
      DO ir=1,NRADCOR
         idgea=IRADCORR(ir)
         ih=0    ! index in the HEP record
 20      ih=ih+1
         IF(ITPTRA(1,ih).EQ.idgea) CALL PHOTOS(ih)  ! PHOTOS may add particles below ih in the list
         IF(ih.LT.NHEP) GO TO 20
      ENDDO
C
C---    Refill local event from PHOTOS output
C
      NTRA = NHEP
      DO i=1,NHEP
         ITPTRA(2,i) = ISTHEP(i)
         ITPTRA(3,i) = IDHEP(i)
         ITPTRA(4,i) = JMOHEP(1,i)
         ITPTRA(5,i) = JDAHEP(1,i)
         ITPTRA(6,i) = JDAHEP(2,i)
         DO j=1,3
            PTRA(j,i) = PHEP(j,i)
         ENDDO
         AMTRA(i) = PHEP(5,i)
C
C---        Restore the IDG from KF
C
         DO j=1,nrc
            IF(ITPTRA(3,i).EQ.kfrc(j)) THEN
               ITPTRA(1,i)=igrc(j)
               GO TO 100
            ENDIF
         ENDDO
 100     CONTINUE
      ENDDO
C
C      IF(IEVENT.LT.NPRIEV) THEN
C         WRITE(6,*) ' ==== After PHOTOS '
C         CALL BG_PRI(1)
C      ENDIF
C
 999  CONTINUE
C      write(6,*) ebeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C

