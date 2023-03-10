* $Header:$
* $Log:$
*
      SUBROUTINE OMROTS2(VAZ,ROT)
C
C     ******************************************************************
C     *                                                                *
C     *       Fill a rotation matrix VA=ROT*VB                         *
C     *  INPUT: VAZ  - 3 vector in frame A, VAZ_norm=ROT*(0,0,1)       *
C     *                X in the radial direction of VAZ                *
C     *                 If VAZ = Z   unit rotation matrix              *
C     * OUTPUT: ROT  - rotation matrix                                 *
C     *                                                                *
C     *    ==>Called by : kinematics programs                          *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
C      REAL     
      DOUBLE PRECISION
     +    VAZ(3),ROT(3,3)
C
      INTEGER   i,j
C      REAL 
      DOUBLE PRECISION
     +     totz,vzt
     +    ,van(3,3)        !  (j,i) j-th projection of the i-th ort (X',Y',Z')
C
C     ------------------------------------------------------------------
C
      DO i=1,3
         ROT(i,i)=1.
         DO j=1,3
            IF(j.NE.i) ROT(j,i)=0.
         ENDDO
      ENDDO
      totz=SQRT(VAZ(1)**2+VAZ(2)**2+VAZ(3)**2)
      IF(totz.GT.0.) THEN
         DO i=1,3
            van(i,3)=VAZ(i)/totz
         ENDDO
         vzt=SQRT(van(1,3)**2+van(2,3)**2)
         IF(vzt.GT.0.) THEN   ! VAZ is not parallel to Z
C            X' projection on XY is along vzt
            DO i=1,2
               van(i,1)=van(3,3)*van(i,3)/vzt
            ENDDO
            van(3,1)=-vzt
C---          Find Y'
            van(1,2)= van(2,3)*van(3,1)-van(3,3)*van(2,1)
            van(2,2)=-van(1,3)*van(3,1)+van(3,3)*van(1,1)
            van(3,2)= van(1,3)*van(2,1)-van(2,3)*van(1,1)
C            write(6,FMT='(3X,3F10.6)') van
C            CALL OVECPROD(van(1,3),van(1,1),van(1,2))  ! find Y'
C            write(6,FMT='(3X,3F10.6)') van
C
            DO i=1,3
               DO j=1,3
                  ROT(i,j)=van(j,i)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
C
      RETURN
      END
