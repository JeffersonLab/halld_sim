      SUBROUTINE OMDECA2(P0,AM,IHEL,PZFR,POUT)
C
C---   2-body decay
C
C---   Input: P0 - initial 4-vector, P0**2 - mass(energy) of the initial state,
C                  defined in the "LAB" frame
C             AM(1:2) - masses of the products
C             IHEL: =0  - cos(th) (to P0 direction) uniform of the 1-st particle
C                         Decay in the Helicity or Gottfried-Jackson frames 
C             PZFR - momentum (in the same frame as P0) defining the Z axis in the Helicity or GJ frame 
C                    after boosting to the decay frame
C             IHEL  <0 - Z is opposite to PZFR in the decay frame (if the recoil is used for the Helicity frame)        
C         ABS(IHEL) =1  - sin(th)**2=(1-cos(th)**2) - for 1-->0+0     rho
C                   =2  - (1+cos(th)**2)            - for 2-->1/2+1/2 psi
C
C      Output: POUT(1:4,1:2) - the secondary 4-momenta
C              
C

      IMPLICIT NONE
      REAL P0(4),AM(2),PZFR(4),POUT(4,2)
      INTEGER IHEL
C
      REAL RNDM,ORNDPOLY
C
      REAL pp(4,2)  ! CM, Z looks along P0
     +    ,pcm(4,2) ! CM, lab orientation
     +    ,bet(4)
     +    ,pzfcm(4) ! the Helicity/GJ vector in the decays frame
     +    ,ecms,ecm,epf1,ppf,ct,st,phi
     +    ,rot(3,3),poly(10),xlim(2),pzm,twopi,vtmp(3),pzf(4)
      INTEGER i,j
C
C--- 
C
C      write(6,*) P0, PZFR
      DO i=1,2
         DO j=1,4
            POUT(j,i)=0.
         ENDDO
      ENDDO
C
      ecms=P0(4)**2-P0(1)**2-P0(2)**2-P0(3)**2
      IF(ecms.LE.0.) THEN
         WRITE(6,*) ' *** OMDECA2 space-like initial vector ',ecms,P0
         GO TO 999
      ENDIF
      ecm=SQRT(ecms)
      IF(ecm.LE.AM(1)+AM(2)) THEN
         WRITE(6,*) ' *** OMDECA2 below threshold ',ecm,AM
         GO TO 999
      ENDIF
C
      epf1=(ecms+AM(1)**2-AM(2)**2)/2./ecm
      ppf=SQRT(epf1**2-AM(1)**2)
C
C---       Boost to CM
C
      bet(4)=p0(4)/ecm
      DO j=1,3
         bet(j)=P0(j)/P0(4)
      ENDDO
C
      IF(IHEL.EQ.0) THEN
         ct=-1.+2.*RNDM(ct)
      ELSE
         xlim(1)=-1.
         xlim(2)= 1.
         IF(ABS(IHEL).EQ.1) THEN
            poly(1)= 1.         ! 1-ct**2    rho  1 --> 0 0
            poly(2)= 0.         ! 
            poly(3)=-1.         ! 
         ELSE IF(ABS(IHEL).EQ.2) THEN
            poly(1)= 1.         ! 1+ct**2    jpsi 1 --> 1/2 1/2
            poly(2)= 0.         ! 
            poly(3)= 1.         ! 
         ENDIF
         ct=ORNDPOLY(2,poly,xlim)         
         IF(ct.LT.-2.) THEN
            WRITE(6,*) ' *** OMDECA2 ct= ',ct,ecm,IHEL
            GO TO 999
         ENDIF
         IF(ABS(ct).GT.1.) THEN
            WRITE(6,*) ' *** OMDECA2 err ct= ',ct,ecm,IHEL
         ENDIF
      ENDIF
C
      twopi=ACOS(0.)*4.
C
      st=SQRT(1.-ct**2)
      phi=twopi*RNDM(st)
C      write(6,*) ' ecm,ct,phi,P0 ',ecm,ct,phi,P0
C
C---    2-body
C
      pp(4,1)=epf1
C     
      pp(1,1)=ppf*st*COS(phi)
      pp(2,1)=ppf*st*SIN(phi)
      pp(3,1)=ppf*ct
C     
      DO j=1,3
         pp(j,2)=-pp(j,1)
      ENDDO
      pp(4,2)=ecm-pp(4,1)
C
C---    For isotropic decays - no rotation
C
      DO i=1,3
         rot(i,i)=1.
         DO j=1,3
            IF(j.NE.i) rot(j,i)=0.
         ENDDO
      ENDDO
C
      IF(IHEL.NE.0) THEN
C
C---      Rotate to the frame where Z goes along PZFR in the decay frame
C
         DO i=1,4
            pzf(i)=PZFR(i)
         ENDDO
         IF(pzf(4).LE.0.) THEN
            WRITE(6,*) ' *** ERROR in OMDECA2 IHEL PZFR =',IHEL,PZFR
            GO TO 100
         ENDIF
C---       Boost PZFR to the decay frame
         CALL GLOREN(bet,pzf(1),pzfcm(1))
C
         pzm=SQRT(pzfcm(1)**2+pzfcm(2)**2+pzfcm(3)**2)
C         write(6,*) ' pzfr ',pzfr,pzm
         IF(pzm.GT.0.00001) THEN
C
            DO i=1,3
               vtmp(i)=0.
               IF(IHEL.LT.0.) pzfcm(i)=-pzfcm(i)
            ENDDO
            CALL OMROTS1(pzfcm,vtmp,rot)
C            CALL OMROTS(P0,rot) ! old - wrong frame based on P0
C            write(6,1111) (pzfcm(i)/pzm,(rot(j,i),j=1,3),i=1,3) ! !!!
C 1111       FORMAT(2X,E12.4,5X,3E12.4) ! !!!
         ELSE
            WRITE(6,*) ' *** ERROR in OMDECA2 PZFR boosted ',pzfcm
            GO TO 100
         ENDIF
C         write(6,FMT='(A4,4F10.4)') 'p0= ' ,(P0 (j),j=1,4) 
C         write(6,FMT='(3F10.4)') ct,pp(3,1)/ppf
C     +        ,(pcm(1,1)*P0(1)+pcm(2,1)*P0(2)+pcm(3,1)*P0(3))/ppf/pzm 
      ENDIF
C
 100  CONTINUE
C
C---       Rotate and boost to the original frame
C
      DO j=1,3
         bet(j)=-bet(j)
      ENDDO
      DO i=1,2
         CALL OMROTV1(0,pp(1,i),rot,pcm(1,i))
         pcm(4,i)=pp(4,i)
         CALL GLOREN(bet,pcm(1,i),POUT(1,i))
      ENDDO
C
C
C      DO j=1,3
C         bet(j)=-bet(j)
C      ENDDO
C      CALL GLOREN(bet,POUT(1,1),pcm(1,1))
C      ct=(pcm(1,1)*P0(1)+pcm(2,1)*P0(2)+pcm(3,1)*P0(3))
C     +  /SQRT(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2)
C     +  /SQRT(P0(1)**2+P0(2)**2+P0(3)**2)
C      write(6,*) ' ct= ',ct
C
C      write(6,*) 'p0=',p0,ecm
C      write(6,FMT='(3F10.4)') rot 
C      write(6,FMT='(4F10.5)') bet 
C      write(6,FMT='(A4,4F10.4)') ('pp= ' ,(pp (j,i),j=1,4),i=1,2) 
C      write(6,FMT='(A4,4F10.4)') ('pcm=' ,(pcm(j,i),j=1,4),i=1,2) 
C      write(6,FMT='(A4,4F10.4)')'dif=',(P0(j)-POUT(j,1)-POUT(j,2),j=1,4)
C      write(6,FMT='(A4,4F10.4)') ('lab=' ,(POUT(j,i),j=1,4),i=1,2) 
C
 999  RETURN
      END

