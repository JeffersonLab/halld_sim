      REAL FUNCTION GET_STRANGE(IFL)
C
C--       Count events with at least one strange particle in the final state 
C
      IMPLICIT NONE
      INTEGER IFL
C
      INCLUDE ?
C
      INTEGER ip,j,ifound
      INTEGER mxstr
      PARAMETER (mxstr=10)
      INTEGER kfs(mxstr)
      DATA kfs/130,321,310,3122,3222,3212,3112,3322,3312,3334/
C
      ifound=0
C
      DO ip=1,NP
         DO j=1,mxstr
            IF(kfs(j).EQ.ITYP(3,ip)) ifound=ifound+1
         ENDDO
      ENDDO
C
      GET_STRANGE=ifound
C
      RETURN
      END
