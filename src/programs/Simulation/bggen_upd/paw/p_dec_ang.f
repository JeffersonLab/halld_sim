      REAL FUNCTION P_DEC_ANG(IFL,K1,K2,K3)
C
C---      Kinematic variables:
C
C        IFL=1 - t  (beam -->K1 : (K1-beam)**2)
C           =2 - cos(th) of K1 in CM
C           =3 - cos(th) of K1 in CM of K2 with respect to the K2 direction
C           =4 - cos(th) of K1 in CM of K2 in the helicity frame, direct
C           =5 - cos(th) of K1 in CM of K2 in the helicity frame, through CM
C           =6 - cos(th) of K1 in CM of K2 in the GJ frame, direct
C           =7 - cos(th) of K1 in CM of K2 in the GJ frame, through CM
C           =8 - phi in the particle momentum transverse plane, X along K2
C         K3>0 - recoil (for the helicity frame)
C
      IMPLICIT NONE
      INTEGER IFL,K1,K2,K3,K4
C
      INCLUDE ?
C
      INTEGER i,j,kf1,kf2
C      REAL 
      DOUBLE PRECISION
     +     var,qq,en1,en2,dir(3),p1(5),p2(5),p3(5),pp,px(4),pcm(4)
     +    ,cthdec(3,2) ! (j,i) i=1 - direct K2, =2 - ->CM->K2; j =1 to p1 (wrong), =2 - Helicit, =3 - GJ 
     +    ,ptar(5),pbeam(5),ptmp(5)
     +    ,betf2(5),p1f2(5),dff2(3,2)  ! bet - boost to K2 CM, p1p2 - K1 in cm2, dff2(3,2) - directions of Z for Helicity and GJ  
     +    ,betcm(4),pbcm(5),p1cm(5),p2cm(5),p3cm(5)  ! in CM
     +    ,betcm2(4),pbcm2(5),p1cm2(5),p3cm2(5),dfcm2(3,2)  ! in CM
     +    ,rotf2(3,3),p1f2n(3),p1f2nr(3)
     +    ,phidec,vtmp(3)
C
      P_DEC_ANG=-20.
      IF(K1.LT.1.OR.K1.GT.NP) GO TO 999
      IF(K2.LT.1.OR.K2.GT.NP) GO TO 999
      IF(K3.LT.1.OR.K3.GT.NP) GO TO 999
C
      DO j=1,3
         ptar(j)=PIN(j,2)
      ENDDO
      ptar(5)=SQRT(ptar(1)**2+ptar(2)**2+ptar(3)**2)
      ptar(4)=SQRT(ptar(5)**2+AMIN(2)**2)
      DO j=1,3
         pbeam(j)=PIN(j,1)
      ENDDO
      pbeam(5)=SQRT(pbeam(1)**2+pbeam(2)**2+pbeam(3)**2)
      pbeam(4)=SQRT(pbeam(5)**2+AMIN(1)**2)
C
      DO j=1,4
         pcm(j)=pbeam(j)+ptar(j)
      ENDDO
      DO j=1,3
         betcm(j)=pcm(j)/pcm(4)
      ENDDO
      betcm(4)=1./SQRT(1.-betcm(1)**2-betcm(2)**2-betcm(3)**2)
C      write(6,*) 'betcm ', betcm
C
      qq=0.
      DO j=1,3
         p1(j)=POUT(j,K1)
         qq=qq+p1(j)**2
      ENDDO
      p1(4)=SQRT(qq+AM(K1)**2)
      p1(5)=SQRT(qq)
C
      qq=0.
      DO j=1,3
         p2(j)=POUT(j,K2)
         qq=qq+p2(j)**2
      ENDDO
      p2(4)=SQRT(qq+AM(K2)**2)
      p2(5)=SQRT(qq)
C
      qq=0.
      DO j=1,3
         p3(j)=POUT(j,K3)
         qq=qq+p3(j)**2
      ENDDO
      p3(4)=SQRT(qq+AM(K3)**2)
      p3(5)=SQRT(qq)
C      write(6,*) 'p1 ', p1
C      write(6,*) 'p2 ', p2
C      write(6,*) 'p3 ', p3
C
C---      Direct boost to K2
C
      DO j=1,3
         betf2(j)=p2(j)/p2(4)
      ENDDO
      betf2(4)=1./SQRT(1.-betf2(1)**2-betf2(2)**2-betf2(3)**2)
C      write(6,*) 'betf2 ', betf2
C
      CALL GLORENF(betf2(1),p1(1),p1f2(1))
      CALL SCALPR(3,p1f2(1),p2(1),qq,cthdec(1,1))

      CALL GLORENF(betf2(1),p3(1),ptmp(1))
      CALL SCALPR(3,p1f2(1),ptmp(1),qq,cthdec(2,1))
      cthdec(2,1)=-cthdec(2,1)
      DO j=1,3
         dff2(j,1)=ptmp(j)/ptmp(5)
         p1f2n(j)=p1f2(j)/p1f2(5)
      ENDDO

      CALL GLORENF(betf2(1),pbeam(1),ptmp(1))
      CALL SCALPR(3,p1f2(1),ptmp(1),qq,cthdec(3,1))
      DO j=1,3
         dff2(j,2)=ptmp(j)/ptmp(5)
      ENDDO
C      CALL GLORENF(betf2(1),p2(1),ptmp(1))
C      write(6,*) ' p1 a ',ptmp
C      CALL OMROTS2(p2,rotf2)
      vtmp(1)=p2(1)
      vtmp(2)=p2(2)
      vtmp(3)=0.
      CALL OMROTS1(dff2(1,2),vtmp,rotf2)
      CALL OMROTV1(1,p1f2n,rotf2,p1f2nr)
      qq=SQRT(p1f2nr(1)**2+p1f2nr(2)**2)
      cthdec(1,2)=p1f2nr(3)
      phidec=0.
      IF(qq.GT.0.) phidec=ACOS(p1f2nr(1)/qq)*180./(ACOS(0.)*2.)
      IF(p1f2nr(2).LT.0.) phidec=360.-phidec
C
C---      Boost to CM, then to K2
C
      CALL GLORENF(betcm,pbeam,pbcm)
      CALL GLORENF(betcm,p1,p1cm)
      CALL GLORENF(betcm,p2,p2cm)
      CALL GLORENF(betcm,p3,p3cm)

      DO j=1,3
         betcm2(j)=p2cm(j)/p2cm(4)
      ENDDO
      betcm2(4)=1./SQRT(1.-betcm2(1)**2-betcm2(2)**2-betcm2(3)**2)
      CALL GLORENF(betcm2,p1cm,p1cm2)
      CALL GLORENF(betcm2,p3cm,p3cm2)
      CALL SCALPR(3,p1cm2(1),p3cm2(1),qq,cthdec(2,2))
      cthdec(2,2)=-cthdec(2,2)
      DO j=1,3
         dfcm2(j,1)=p3cm2(j)/p3cm2(5)
      ENDDO
C
      CALL GLORENF(betcm2,pbcm,pbcm2)
      CALL SCALPR(3,p1cm2(1),pbcm2(1),qq,cthdec(3,2))
      DO j=1,3
         dfcm2(j,2)=pbcm2(j)/pbcm2(5)
      ENDDO
C      CALL GLORENF(betcm2(1),p2cm(1),ptmp(1))
C      write(6,*) ' p1 b ',ptmp
C
C      write(6,1000) cthdec
C 1000 format(3X,3F12.5)
C      write(6,FMT='(2X,A10,3X,F10.3)') 'phidec',phidec
C      write(6,*) ' -------------'
C      
      var=-20.
C
      IF(IFL.EQ.1) THEN
C
         var=AM(K1)**2+AMIN(1)**2-2.*p1(4)*pbeam(4)
         DO j=1,3
            var=var+2.*p1(j)*pbeam(j)
         ENDDO
C
      ELSE IF(IFL.EQ.2) THEN
         CALL SCALPR(3,p1cm,pbcm,qq,var)
      ELSE IF(IFL.EQ.3) THEN
         var=cthdec(1,1)
      ELSE IF(IFL.EQ.4) THEN
         var=cthdec(2,1)
      ELSE IF(IFL.EQ.5) THEN
         var=cthdec(2,2)
      ELSE IF(IFL.EQ.6) THEN
         var=cthdec(3,1)
      ELSE IF(IFL.EQ.7) THEN
         var=cthdec(3,2)
      ELSE IF(IFL.EQ.8) THEN
         var=phidec
      ELSE IF(IFL.EQ.9) THEN
         var=p1f2nr(1)
      ELSE IF(IFL.EQ.10) THEN
         var=p1f2nr(2)
      ENDIF

C
      P_DEC_ANG=var
C
 999  RETURN
C
      END
C
      SUBROUTINE SCALPR(N,A,B,SP,SPN)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION  A(*),B(*),SP,SPN
      DOUBLE PRECISION  qa,qb
      INTEGER i
C
      SP=0.
      qa=0.
      qb=0.
      DO i=1,N
         SP=SP+A(i)*B(i)
         qa=qa+A(i)**2
         qb=qb+B(i)**2
      ENDDO
      qa=SQRT(qa)
      qb=SQRT(qb)
      SPN=0.
      IF(qa.GT.0.AND.qb.GT.0) SPN=SP/qa/qb
C
      RETURN
      END
C
C      INCLUDE 'efmass.f'
      INCLUDE 'glorenf.f'
      INCLUDE 'omrots2.f'
      INCLUDE 'omrotv1.f'
      INCLUDE 'omrots1.f'
      INCLUDE 'ovecprod.f'
