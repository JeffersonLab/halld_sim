* $Header:$
* $Log:$
*
      SUBROUTINE OMROTV1(IFL,V1,ROT,V2)
C
C     ******************************************************************
C     *                                                                *
C     *       Vector rotation V1 ==> V2 using ROT matrix               *
C     *                                                                *
C     *    ==>Called by : OMKINE                                       *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER  IFL     ! =0 - forward v2=rot*v1
C                        >0 - backward v2=rot^(-1)*v1  ( rot^(-1)=rot^T
      REAL     V1(3),ROT(3,3),V2(3)
C
      INTEGER   i,j
      REAL r
C
C     ------------------------------------------------------------------
C
      DO i=1,3
         V2(i)=0.
         DO j=1,3
            IF(IFL.EQ.0) THEN
               r=ROT(j,i)
            ELSE
               r=ROT(i,j)
            ENDIF
            V2(i)=V2(i)+r*V1(j)
         ENDDO
      ENDDO
C
      RETURN
      END






