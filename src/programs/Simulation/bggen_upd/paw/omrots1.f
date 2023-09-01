* $Header:$
* $Log:$
*
      SUBROUTINE OMROTS1(VAZ,VAX0,ROT)
C
C     ******************************************************************
C     *                                                                *
C     *       Fill a rotation matrix VA=ROT*VB                         *
C     *  INPUT: VAZ  - 3 vector in frame A, VAZ=ROT*(0,0,Z) (Z>0)      *
C     *         VAX0 - 3 vector in frame A, VAX=ROT*(X,0,0) (X>0)      *
C     *               if VAX0**2=0 - select X in B closer to           *
C     *                             the smallest projection of VAZ     *
C     *               if VAX0**2>0 - select Y in B as VAZxVAX0         *
C     * OUTPUT: ROT  - rotation matrix                                 *
C     *                                                                *
C     *    ==>Called by : kinematics programs                          *
C     *                                                                *
C     ******************************************************************
C
      IMPLICIT NONE
C      REAL     
      DOUBLE PRECISION
     +       VAZ(3),VAX0(3),ROT(3,3)
C
      INTEGER   i,j,m
     +         ,ky        ! >0 - Y defined using VAX
     +         ,kord(3)   ! order of V projections (absolute values) in increasing order 
C      REAL
      DOUBLE PRECISION 
     +     totz,totx,q
     +    ,vazn(3),vazna(3)   ! normalized VAZ
     +    ,vax0n(3)        !    normalized VAX0
     +    ,vaxn(3)         ! vaxn=R*(1,0,0) - projections of XB on A   
     +    ,vayn(3)         ! vayn=R*(0,1,0) - projections of YB on A   
C
C     ------------------------------------------------------------------
C
      totz=SQRT(VAZ(1)**2+VAZ(2)**2+VAZ(3)**2)
      IF(totz.GT.0.) THEN

         DO i=1,3
            DO j=1,3
               ROT(j,i)=0.
            ENDDO
         ENDDO
C
         ky=0
         totx=SQRT(VAX0(1)**2+VAX0(2)**2+VAX0(3)**2)
         IF(totx.GT.0.) ky=1 
         DO i=1,3
            vazn(i)=VAZ(i)/totz
            vazna(i)=ABS(vazn(i))
            IF(ky.GT.0) vax0n(i)=VAX0(i)/totx
         ENDDO
C
C---     Define X-Y vectors normal to vazn
C
         IF(ky.GT.0) THEN
C---       Use VAX0 to define Y 
            CALL OVECPROD(vazn,vax0n,vayn)
            q=SQRT(vayn(1)**2+vayn(2)**2+vayn(3)**2)
            IF(q.GT.1.E-10) THEN
               DO i=1,3
                  vayn(i)=vayn(i)/q
               ENDDO
               CALL OVECPROD(vayn,vazn,vaxn)
            ELSE
               WRITE(6,1000) VAZ,VAX0,vayn
 1000          FORMAT(' *** OMROTS1 error: VAZ and VAX are parallel '
     +                ,3(3E12.4,4X))
               ky=0
            ENDIF
         ENDIF        
C
         IF(ky.EQ.0) THEN    !  VAX0 is not defined
            DO i=1,3
               vayn(i)=0.
            ENDDO
            CALL SORTIN(3,vazna(1),kord(1))
            q=vazn(kord(1))**2+vazn(kord(2))**2
C            write(6,*) kord,q
            IF(q.LT.1.E-10) THEN
C---           VAZ is along one of the axis
               m=kord(3)-1         !  Y direction
               IF(m.LT.1) m=m+3
            ELSE
C---           VAZ is not along any of the axis 
               m=kord(1)           ! along the smallest VAZ component in A - Y in B
            ENDIF 
            vayn(m)=1. 
            CALL OVECPROD(vayn,vazn,vaxn)
C            write(6,*) vayn,'   ',vaxn
            q=SQRT(vaxn(1)**2+vaxn(2)**2+vaxn(3)**2)
            DO i=1,3
               vaxn(i)=vaxn(i)/q
            ENDDO
            CALL OVECPROD(vazn,vaxn,vayn)
         ENDIF
         DO i=1,3
            ROT(1,i)=vaxn(i)       ! R*(1,0,0)=vaxn
            ROT(2,i)=vayn(i)       ! R*(0,1,0)=vayn
            ROT(3,i)=vazn(i)       ! R*(0,0,1)=vazn
         ENDDO
C
      ELSE
         DO i=1,3
            DO j=1,3
               IF(j.EQ.i) THEN
                  ROT(j,i)=1.
               ELSE
                  ROT(j,i)=0.
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
