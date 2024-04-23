      SUBROUTINE PYTH_PRI(IFL,AMLIM)
C
C---   Prints all PYTHIA particles with mass < AMLIM
C---   IFL=-1 - no dacays
C---       -2 - with decays
C
      IMPLICIT NONE
      INTEGER IFL
      REAL AMLIM
C      
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C
      DOUBLE PRECISION PMAS,PARF,VCKM,BRAT
      INTEGER          KCHG,MDCY,MDME,KFDP
C
      INTEGER kc,kf,i,ifd,kdc,ndc,idc,kfd,lenn
      REAL am
      CHARACTER chap*16,chan*16,chad*16
C
      IF(IFL.GT.0) THEN
         WRITE(6,1010)
 1010    FORMAT(3X,'kc',8X,'kf',2X,'dec',2X,'particle',9X
     +         ,'antiparticle',7X,'mass',5X,'width')
         DO kc=1,500
            kf=KCHG(kc,4)
            IF(kf.NE.0) THEN
               CALL PYNAME(kf,chap)
               chan=' '
               IF(MDCY(kc,3).NE.0) CALL PYNAME(-kf,chan)
               IF(kf.GT.100.OR.chap.NE.' ') THEN
                  am=PMAS(kc,1)
C   remove high mass and zero mass except neutrinos and photon
                  IF(am.LT.AMLIM.AND.
     +              (am.GT.0..OR.kf.EQ.12.OR.kf.EQ.14.OR.kf.EQ.22)) THEN 
                     WRITE(6,1020) kc,kf,MDCY(kc,1),chap,chan
     +                 ,(PMAS(kc,i),i=1,2)
 1020                FORMAT(1X,80('='),/,I5,2X,I8,2X,I2,3X,A16,1X,A16
     +                     ,1X,2F9.5)
                     IF(IFL.GE.2) THEN
                        ifd=MDCY(kc,1) ! decays allowed
                        kdc=MDCY(kc,2) ! address
                        ndc=MDCY(kc,3) ! number of decays
                        IF(ndc.EQ.0.OR.kdc.EQ.0) GO TO 100
                        IF(IFL.EQ.2.and.ifd.EQ.0) GO TO 100
                        WRITE(6,FMT='(1X,80(''-''))')
                        DO idc=kdc,kdc+ndc-1
                           WRITE(6,1050,advance='no')
     +                           (MDME(idc,i),i=1,2),BRAT(idc)
 1050                      FORMAT(2X,I3,2X,I3,2X,F9.6,1X)
                           DO i=1,5
                              kfd=KFDP(idc,i)
                              IF(kfd.NE.0) THEN
                                 CALL PYNAME(kfd,chad)
                                 lenn=LEN_TRIM(chad)
                                 WRITE(6,1060,advance='no') kfd
     +                                ,chad(1:lenn) 
 1060                            FORMAT(1X,I7,1X,A)
                              ENDIF
                           ENDDO
                           WRITE(6,FMT='(1X)')
                        ENDDO
 100                    CONTINUE
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         WRITE(6,FMT='(1X,80(''-''))')
C
      ENDIF
C
      RETURN
      END
