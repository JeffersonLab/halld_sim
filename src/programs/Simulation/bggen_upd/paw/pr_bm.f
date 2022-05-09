      SUBROUTINE PR_BM(E1,E2,NB)
      IMPLICIT NONE 
C---     Print a "beam spectrum"
C
      REAL E1,E2
      INTEGER NB
C
      INTEGER i,np
      REAL e,de
      REAL sp(5000),ep(5000)
C
      de=(E2-E1)/NB
      np=0
      DO i=1,NB
         e=E1+(i-0.5)*de
         IF(e.GT.(E1+E2)/2..OR.MOD(i-1,4).EQ.0) THEN
            np=np+1
            ep(np)=e
            sp(np)=E2/e
         ENDIF
      ENDDO
      WRITE(6,1000) (sp(i),i=1,np)
 1000 FORMAT(8X,10F8.4)
      WRITE(6,*) '  np=',np
      WRITE(6,1000) (ep(i),i=1,np)
C
      RETURN
      END
