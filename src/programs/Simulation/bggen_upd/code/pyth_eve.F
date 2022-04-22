      SUBROUTINE PYTH_EVE(IERR)
C
C---  Simulates 1 PYTHIA event 
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
C
      INTEGER PYK,PYCOMP
      DOUBLE PRECISION PYP
      EXTERNAL PYK,PYP,PYCOMP
C
      REAL      beamen
C
      INTEGER i,j,nlnd,ilnd,kf,kc,ks
C
C     ------------------------------------------------------------------
C
C---      Variable energy?
C
      beamen=PIN(3,1)
      IF(NHBEA.GT.1) THEN
         CALL PARP_INI(171,beamen/EPH_LIM(2)) ! the relative energy for this event
      ENDIF
C
      CALL PYEVNT
      CALL PYEDIT(15)  ! filter out some intermediate entries
      nlnd=PYK(0,1)
      NTRA=nlnd
C
      DO ilnd=1,MIN(nlnd,MXTRA)
         DO i=1,5
            ITPTRA(i+1,ilnd)=PYK(ilnd,i)
         ENDDO
         ks=ITPTRA(2,ilnd)
         kf=ITPTRA(3,ilnd)
         kc=PYCOMP(kf)
         IF(kf.LT.0) kc=-kc
C         write(6,*) ilnd,kf,kc,MXPKC,KCGEAN(kc)
         IF(ABS(kc).LE.MXPKC) THEN
            ITPTRA(1,ilnd)=IDG_KC(kc)
            MULT_KC(kc)=MULT_KC(kc)+1    ! multiplicity of this particle type
            IF(ks.LT.11) NODC_KC(kc)=NODC_KC(kc)+1    ! number of "undecayed" particles of this type
            IF(MULT_KC(kc).EQ.1) THEN
C               CALL PYNAME(kf,CH_KC(kc)) ! first occurrence - initialize the particle name
               KF_KC(kc)=kf
            ENDIF
         ELSE
            ITPTRA(1,ilnd)=0
         ENDIF
         DO i=1,3
            PTRA(i,ilnd)=REAL(PYP(ilnd,i))
         ENDDO
         AMTRA(ilnd)=REAL(PYP(ilnd,5))
C
      ENDDO
C
      IEVPROC=0
      IERR=0
 999  CONTINUE
C
      END

