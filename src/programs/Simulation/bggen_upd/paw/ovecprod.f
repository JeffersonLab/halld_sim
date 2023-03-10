      SUBROUTINE OVECPROD(V1,V2,VR)
C
C---   Vector product of 3-dim vectors VR=V1xV2
C
      IMPLICIT NONE
C      REAL 
      DOUBLE PRECISION
     +        V1(3),V2(3),VR(3)
C
      VR(1)= V1(2)*V2(3)-V1(3)*V2(2)
      VR(2)=-V1(1)*V2(3)+V1(3)*V2(1)
      VR(3)= V1(1)*V2(2)-V1(2)*V2(1)
C
      RETURN
      END

