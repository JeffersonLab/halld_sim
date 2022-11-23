* $Header: /afs/cern.ch/exp/compass/src/cvs/comgeant/code/src/omgbatch/ompro/sortin.F,v 1.1.1.1 1997/06/02 17:39:56 fdr Exp $
* $Log: sortin.F,v $
* Revision 1.1.1.1  1997/06/02 17:39:56  fdr
* Comgeant Monte-Carlo
*
* Revision 3.2.0.1  1996/11/07 19:23:30  las
* First CVS version.
*
*CMZ :          15/12/91  22.19.32  by  E.Chudakov
*-- Author :    E.Chudakov   07/02/91
C
C======================================================================
C
      SUBROUTINE SORTIN(N,X,IPOIN)
      IMPLICIT NONE
      INTEGER N,IPOIN(*)
      REAL X(*)
C
C     Sorting of N real values X. IPOIN - numbers of values X in
C      increasing order.
C
      INTEGER i,j1,k,j2
      REAL xm
C
      IF(N.LE.0) RETURN
C
      DO i=1,N
         IPOIN(i)=i
      ENDDO
C
      IF(N.EQ.1) RETURN
C
      DO i=1,N-1
         j1=IPOIN(i)
         xm=X(j1)
         DO k=i+1,N
            j2=IPOIN(k)
            IF(X(j2).LT.xm) THEN
               IPOIN(i)=j2
               IPOIN(k)=j1
               xm=X(j2)
               j1=j2
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
