      SUBROUTINE BG_PRI(IPRI)
C
C---  Prints the event record
C
      IMPLICIT NONE
      INTEGER IPRI
C
      INCLUDE 'bg_evec.inc'
      INCLUDE 'bg_partc.inc'
C
      INTEGER i,j,lout,it,ig,ip,kf,kc
      CHARACTER cent(2)*6,cproc*16,cpnam*16
C
C     ------------------------------------------------------------------
C
      lout=6
      cent(1)='beam  '
      cent(2)='target'
C
      IF(IPRI.NE.0) THEN
         cproc=CPROCNAM
         WRITE(lout,1000) IEVENT,IEVPROC,CPROCNAM
 1000    FORMAT(1X,80('='),/,' Event ',I6,'  Process=',I4,3X,A)
         WRITE(lout,1002)
 1002    FORMAT(11X,'ID',8X,'KF',23X,'mass',6X
     +        ,'Px',6X,'Py',6X,'Pz',6X,'name')
         DO i=1,2
            cpnam='?'
            ig=ITPIN(1,i)
            kf=ITPIN(2,i)
            IF(ig.GT.0) cpnam=CH_PART(ig)
            WRITE(lout,1005) 
     +       cent(i),ig,kf,AMIN(i),(PIN(j,i),j=1,3),cpnam
 1005       FORMAT(1X,A6,3X,I3,5X,I5,19X,F8.4,3X,3F8.3,3X,A)
         ENDDO
         IF(NTRA.GT.0) THEN
            WRITE(lout,1006) 
 1006       FORMAT(1X,80('-'),/
     +            ,'   #',3X,' ID',4X,'KS',5X,'KF',3X,'orig','  decays'
     +            ,8X,'mass',6X,'Px',6X,'Py',6X,'Pz',6X,'name')   
         ENDIF
         DO it=1,NTRA
            cpnam='?'
            kf=ITPTRA(3,it)
            IF(kf.NE.0) THEN
               DO i=1,NPADEF
                  ig=IDGDEF(i)
                  IF(ig.GT.0.AND.ig.LE.MXPART) THEN
                     IF(KF_PART(ig).EQ.kf) cpnam=CH_PART(ig)
                  ENDIF
               ENDDO
               IF(cpnam.EQ.'?') THEN ! try PYTHIA particle names
                  CALL PYNAME(kf,cpnam)
C                  DO kc=-MXPKC,MXPKC
C                     IF(kc.NE.0.AND.KF_KC(kc).EQ.kf) THEN
C                        cpnam=CH_KC(kc) 
C                     ENDIF
C                  ENDDO
               ENDIF
            ENDIF
            WRITE(lout,1010) 
     +      it,(ITPTRA(j,it),j=1,6),AMTRA(it),(PTRA(j,it),j=1,3),cpnam
 1010       FORMAT(1X,I3,3X,I3,I6,I7,3X,3I4,4X,F8.4,3X,3F8.3,3X,A)
         ENDDO
      ENDIF
C
      END
C

