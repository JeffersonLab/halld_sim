      SUBROUTINE PYTH_GETP(KF,KC,CNAM,AMASS,WIDTH,NDEC,IFDEC)
C
C---   Get the main parameters of the particle KF
C---   INPUT:   KF
C---   OUTPUT:  KC <= 0 - no such partice (intrinsic code)
C               CNAM - name
C               AMASS
C               WIDTH
C               NDEC - number of decays
C               IFDEC >0 decays allowed in PYTHIA
C               IERR >0 - no such particle
C
      IMPLICIT NONE
      INTEGER KF,KC,NDEC,IFDEC,IERR
      REAL AMASS,WIDTH
      CHARACTER CNAM*(*)
C      
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C
      DOUBLE PRECISION PMAS,PARF,VCKM,BRAT
      INTEGER          KCHG,MDCY,MDME,KFDP
C
      INTEGER  PYCOMP
      EXTERNAL PYCOMP
C
      KC=PYCOMP(KF)
      IF(KC.LE.0) GO TO 999
      CALL PYNAME(KF,CNAM)
      AMASS=PMAS(KC,1)
      WIDTH=PMAS(KC,2)
      NDEC =MDCY(kc,3)
      IFDEC=MDCY(kc,1)
C
 999  RETURN
      END
