c*********************************************************************************************
c   Four-vectors are of the form (px,py,pz,m,E) - this is a change from the original version 
c   of this program, where the 5th entry was not defined.
c
c	Changes:
c	11-1-91	R.A.Sch. Unbug in ROTATE for rotations about zero-length axes.
c	1-25-94 R.A.Sch. Unbug definition of DIR3 in KIN12
c       8-18-97   "      UNIX version
c      11-12-13   "      Make KIN1R give c.m. momentum back to calling program
c      12-16-14   "      Change BOOST to return the total energy of the boosted 4-vector
c       9-30-15   "      Tweak code to make gfortran on Mac stop complaining
c       6-28-17   "      Make more of the 4-vectors of the form (px,py,pz,mass,energy)
c
      SUBROUTINE EPCM(P1,P2,PCM)
c calculates total momentum and cm energy (invariant mass) of two particles
      REAL P1(5),P2(5),PCM(5)
      AP1S = P1(1)**2+P1(2)**2+P1(3)**2
      AP2S = P2(1)**2+P2(2)**2+P2(3)**2
      APCMS=0.
      DO 100 I=1,3
         PCM(I) = P1(I)+P2(I)
         APCMS  = APCMS+PCM(I)**2
100   CONTINUE
      E=SQRT(AP1S+P1(4)**2)+SQRT(AP2S+P2(4)**2)
      PCM(4) = SQRT(E*E-APCMS)  ! W
      PCM(5) = E
      RETURN
      END

      FUNCTION ABSP(P)
      REAL P(5)
      ABSP=SQRT(P(1)**2+P(2)**2+P(3)**2)
      RETURN
      END

      FUNCTION PSQ(P)
      REAL P(5)
      PSQ=P(1)**2+P(2)**2+P(3)**2
      RETURN
      END

      FUNCTION ETOT(P)
      REAL P(5)
      ETOT = SQRT(P(1)**2+P(2)**2+P(3)**2+P(4)**2)
      if(abs(etot - P(5)).gt.1.e-6)then
         write(6,*)'ETOT routine is not happy',P,ETOT
         call exit()
      endif
      RETURN
      END

      FUNCTION DOT(P1,P2)
      REAL P1(3),P2(3)
      DOT=P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3)
      RETURN
      END

      SUBROUTINE KIN2R(P1,P2,P3,P4,CTHP3,IERR)
c 2 body kinematics S-wave P1+P2 to P3+P4
      REAL P0(5),P1(5),P2(5),P3(5),P4(5)
      CALL EPCM(P1,P2,P0)
      CALL KIN1R(P0,P3,P4,CTHP3,PCM,IERR)
      RETURN
      END

      SUBROUTINE KIN1R(P0,P1,P2,CTHCM,PCM,IERR)
c     Decay kinematics. Random angles
c     P0 - 4-vector (px,py,pz,m,E) of parent particle
c     P1 - first  decay product's 4-vector
c     P2 - second decay product's 4-vector
c     CTHCM - randomly-selected cosine of decay angle with respect to P0 direction
c     PCM - magnitude of the rest-frame P0 momentum of the decay products
c     IERR - signals that an error was detected
c
      REAL P0(5),P1(5),P2(5)
      REAL P0TMP(5),P1TMP(5),P2TMP(5),DTMP(5),BETA(3)
      COMMON/monte/ISEED
c      DATA ISEED /123456/
c     Pick angles relative to velocity of P0
      CTHCM=1.-2.*RAN(ISEED)
      PHI=2.*3.14159*RAN(ISEED)      
c      cthcm = 0.0      !Test hacks
c      phi = 3.1415/2.0
c      phi = 0
      E0=ETOT(P0)
c     velocity of lab frame relative to CM
      do I=1,3
         P0TMP(I) = 0.
         BETA(I)  = -P0(I)/E0
      enddo
      P0TMP(4)=P0(4)
      TMP=1.-CTHCM**2
      IF(TMP.GT.0.) THEN
         STHCM=SQRT(TMP)
      ELSE
         STHCM=0.
      ENDIF
      DTMP(1)=STHCM*COS(PHI)
      DTMP(2)=STHCM*SIN(PHI)
      DTMP(3)=CTHCM
      P1TMP(4)=P1(4)
      P2TMP(4)=P2(4)
c     Let P0 decay at rest with Particle 2 in direction DTMP
      CALL KIN12(P0TMP,P1TMP,P2TMP,DTMP,IERR)
      CALL BOOST(P1TMP,BETA,P1)
      CALL BOOST(P2TMP,BETA,P2)
      PCM = ABSP(P1TMP)
c      write(6,*)'p0  ',p0
c      write(6,*)'p1  ',p1
c      write(6,*)'p2  ',p2
c      write(6,*)'beta',beta
c      write(6,*)'dtmp',dtmp
c      write(6,*)'p1tmp',p1tmp
c      write(6,*)'p2tmp',p2tmp
c      write(6,*)'pcm  ',pcm
      RETURN
      END

      SUBROUTINE KIN1RM(P0,P3,P1,P2,P4,P5,CTHCM,IERR)
c     decay kinematics. Random angles.
c     lets two particles decay in exactly the same way at rest.
c     also returns cosine of P1 and P4 relative to P0.
c     P0 goes to P1 and P2, P3 goes to P4 and P5.
c
      REAL P0(5),P1(5),P2(5),P3(5),P4(5),P5(5)
      REAL P0TMP(5),P1TMP(5),P2TMP(5),P3TMP(5),P4TMP(5)
     1    ,P5TMP(5),DTMP(5)
      REAL BETA0(3),BETA3(3)
      COMMON /monte/ ISEED
c     pick angles relative to velocity of P0 and P3
      CTHCM=1.-2.*RAN(ISEED)
      PHI=2.*3.14159*RAN(ISEED)
      E0=ETOT(P0)
      E3=ETOT(P3)
c     velocity of lab frame relative to CM of P0
      DO 100 I=1,3
      P0TMP(I)=0.
100   BETA0(I)=-P0(I)/E0
      P0TMP(4)=P0(4)
c     velocity of lab frame relative to CM of P3
      DO 200 I=1,3
      P3TMP(I)=0.
200   BETA3(I)=-P3(I)/E3
      P3TMP(4)=P3(4)
      TMP=1.-CTHCM**2
      IF(TMP.GT.0.) THEN
        STHCM=SQRT(TMP)
      ELSE
        STHCM=0.
      ENDIF
      DTMP(1)=STHCM*COS(PHI)
      DTMP(2)=STHCM*SIN(PHI)
      DTMP(3)=CTHCM
      P1TMP(4)=P1(4)
      P2TMP(4)=P2(4)
      P4TMP(4)=P4(4)
      P5TMP(4)=P5(4)
c     let P0TMP and P3TMP decay at rest with particles P2TMP and P4TMP in direction
c     DTMP respectively.
      CALL KIN12(P0TMP,P1TMP,P2TMP,DTMP,IERR)
      CALL KIN12(P3TMP,P4TMP,P5TMP,DTMP,IERR)
      CALL BOOST(P1TMP,BETA0,P1)
      CALL BOOST(P2TMP,BETA0,P2)
      CALL BOOST(P4TMP,BETA3,P4)
      CALL BOOST(P5TMP,BETA3,P5)
      RETURN
      END

      SUBROUTINE KIN22(P1,P2,P3,P4,DIR3,IERR)
c     2 body kinematics P1+P2 goes to P3+P4, P3 in direction DIR3
c     returns momentum in P3 and P4
      REAL P1(5),P2(5),P3(5),P4(5),DIR3(3)
      REAL P0(5)
      CALL EPCM(P1,P2,P0)
      CALL KIN12(P0,P3,P4,DIR3,IERR)
      RETURN
      END

c     Calculates calculates decay of 1 particle to 2 particles
c     P0   = momentum and mass of initial particle
c     DIR3 = directional cosines of particle P3 (not necessarily normalized)
c     P3(4), P4(4) = masses of decay products
c     Outputs 3 vector momenta and total energy in P3 and P4
c
      SUBROUTINE KIN12(P0,P3,P4,DIR3,IERR)
      REAL P0(5),DIR3(3),P3(5),P4(5)
      AP0S=PSQ(P0)
      ADIR=SQRT(dir3(1)*dir3(1) + dir3(2)*dir3(2) + dir3(3)*dir3(3))
c     calculate cosine of angle between P0 and DIR3
      IF(AP0S.NE.0.) THEN
         CALPH=DOT(P0,DIR3)/SQRT(AP0S)/ADIR
      ELSE
         CALPH=1.
      ENDIF
      E0S=P0(4)**2+AP0S
      E0=SQRT(E0S)
      D=E0S+P3(4)**2-P4(4)**2-AP0S
      A=4.*(AP0S*CALPH**2-E0S)
      B=4.*D*SQRT(AP0S)*CALPH
      C=D*D-4.*E0S*P3(4)**2
c     assume highest value P3 if two solutions:
      E=B*B-4.*A*C
      IF(E.GE.0) THEN
         IERR=0
         AP3=(-B-SQRT(E))/(2.*A)
         DO 100 I=1,3
         P3(I)=AP3*DIR3(I)/ADIR
         P4(I)=P0(I)-P3(I)
100      CONTINUE
      ELSE
         IERR=1
      ENDIF
      p3(5) = sqrt(p3(1)**2+p3(2)**2+p3(3)**2+p3(4)**2) !Energy
      p4(5) = sqrt(p4(1)**2+p4(2)**2+p4(3)**2+p4(4)**2) !Energy
      RETURN
      END

c
c     Boosts momentum P1 to frame with velocity BETA3(I)
c     (px,py,pz,m,E) is the form of the '4'-vectors
c 
      SUBROUTINE BOOST(P1,BETA3,P2)
      REAL P1(5),BETA3(3),P2(5),PTMP(5)
c 
      PTMP(4) = P1(4)
      P2(4)   = P1(4)
      BET = sqrt(beta3(1)**2. + beta3(2)**2. + beta3(3)**2.)
      GAM = 1./SQRT(1.-BET**2)
c      write(6,*)'Boost in  ',p1
      E1  = P1(5)
c
c     Rotate to coordinate system with the boost along Z-axis
c     Do the boost, then rotate back...  
c
      CALL ROTATE(P1,BETA3,PTMP)
      PTMP(3) = GAM*(PTMP(3)-BET*E1)
      CALL UNROT(PTMP,BETA3,P2)
      p2(5) = sqrt(p2(1)**2+p2(2)**2+p2(3)**2+p2(4)**2) !Energy
c      write(6,*)'Boost out ',p2
      RETURN
      END

      SUBROUTINE JACOB(P1,BET,P2,THETA,AJAC)
c     Boosts momentum P1 to frame moving with velocity BET along z-axis
c     Also returns jacobian needed for differentical cross section
c     and theta (degrees) angle of P2 w.r.t. Z-axis
      REAL P1(5),P2(5)
      GAM=1./SQRT(1.-BET**2)
      E1=ETOT(P1)
      P2(1)=P1(1)
      P2(2)=P1(2)
      P2(3)=GAM*P1(3)+BET*GAM*E1
      P2(4)=P1(4)
      P2ABS=ABSP(P2)
      BET2=P2ABS/ETOT(P2)
      COS2=P2(3)/P2ABS
      THETA=180./3.1415926536*ACOS(COS2)
c     now compute Jacobian according to CERN notes
      IF(BET.EQ.0) THEN
        AJAC=1.
      ELSE
         GAMSQ=1./(1.-BET*BET)
         AJAC= (1. -BET/BET2 *COS2) /
     >    (GAMSQ*((BET/BET2 -COS2)**2 + (1.-COS2*COS2)/GAMSQ)**1.5)
      ENDIF
      RETURN
      END

      SUBROUTINE ROTATE(P1,DIR3,P2)
c     Rotates vector P1 into direction given by direction cosines DIR3
c     (rotates vector DIR3 into (0,0,1) )
c     Rotate out X-coord, sine and cosine of angle about Z
      REAL P1(3),P2(3),DIR3(3)
      R=SQRT(DIR3(1)**2+DIR3(2)**2)
      IF(R.EQ.0.) THEN
         SINA1=0.
         COSA1=1.
      ELSE
         SINA1=-DIR3(1)/R
         COSA1=DIR3(2)/R
      ENDIF
      YTMP=-SINA1*DIR3(1)+COSA1*DIR3(2)
c     rotate out Y-coord. rotate about X-axis
      R=SQRT(YTMP**2+DIR3(3)**2)
      IF(R.EQ.0.) THEN
         SINA2=0.
         COSA2=1.	!R.A.Sch.        COSA2=0.
      ELSE
         SINA2=-YTMP/R
         COSA2=DIR3(3)/R
      ENDIF
c     now transform
      P2(1)=COSA1*P1(1)+SINA1*P1(2)
      YTMP=-SINA1*P1(1)+COSA1*P1(2)
      P2(2)=COSA2*YTMP+SINA2*P1(3)
      P2(3)=-SINA2*YTMP+COSA2*P1(3)
      RETURN
      END

      SUBROUTINE UNROT(P1,DIR3,P2)
c     Unrotates the rotate of routine ROTATE
c     first work out the two angles
c     rotate out X-coord, sine and cosine of angle about Z
      REAL P1(3),P2(3),DIR3(3)
      R=SQRT(DIR3(1)**2+DIR3(2)**2)
      IF(R.EQ.0.) THEN
         SINA1=0.
         COSA1=1.
      ELSE
         SINA1=-DIR3(1)/R
         COSA1=DIR3(2)/R
      ENDIF
      YTMP=-SINA1*DIR3(1)+COSA1*DIR3(2)
c     rotate out Y-coord. rotate about X-axis
      R=SQRT(YTMP**2+DIR3(3)**2)
      IF(R.EQ.0.) THEN
         SINA2=0.
         COSA2=1.
      ELSE
         SINA2=-YTMP/R
         COSA2=DIR3(3)/R
      ENDIF
c     now transform
      P2(3)=SINA2*P1(2)+COSA2*P1(3)
      YTMP=COSA2*P1(2)-SINA2*P1(3)
      P2(2)=SINA1*P1(1)+COSA1*YTMP
      P2(1)=COSA1*P1(1)-SINA1*YTMP
      RETURN
      END

      SUBROUTINE NEARST(A1,V1,A2,V2,R1,R2,D,IERR)
c
c     find distance and points of closest approach between 2 lines in 3D
c     line one given by R1(i)=A1(i)+t*V1(i)  i=1,3
c     outputs R1,R2 and distance between them, D
c
c     IERR=2 if V1 or V2 are zero
c     IERR=1 if lines are parallel (but "D" is still valid)
      REAL A1(3),V1(3),A2(3),V2(3),R1(3),R2(3)
      V11=V1(1)**2+V1(2)**2+V1(3)**2
      V22=V2(1)**2+V2(2)**2+V2(3)**2
      V12=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
      DET=V12*V12-V11*V22
      IF(V11.LT.1.E-9 .OR. V22.LT.1.E-9) THEN
          IERR=2
          D=1.E10
          GOTO 950
      ENDIF
      ADIF1=A1(1)-A2(1)
      ADIF2=A1(2)-A2(2)
      ADIF3=A1(3)-A2(3)
      C1=ADIF1*V1(1)+ADIF2*V1(2)+ADIF3*V1(3)
      C2=ADIF1*V2(1)+ADIF2*V2(2)+ADIF3*V2(3)
      IF(ABS(DET).LT.1.E-9) THEN
         IERR=1
         T1=0.
         T2=C2/V22
      ELSE
         IERR=0
         T1=(V22*C1-V12*C2)/DET
         T2=(V12*C1-V11*C2)/DET
      ENDIF
      D=0.
      DO 100 I=1,3
         R1(I)=A1(I)+T1*V1(I)
         R2(I)=A2(I)+T2*V2(I)
         D=D+(R1(I)-R2(I))**2
100   CONTINUE
      D=SQRT(D)
950   RETURN
      END

      SUBROUTINE TRAVEL(RIN,DIR,D,ROUT)
c     calculates new position of particle starting at position RIN(i)
c     travelling in direction DIR(i).  Output is ROUT(i)
c     if DIR(i) are directional cosines, input D is distance
c     if DIR(i) is BETA(i), D is time
c     if DIR(i) is un-normalized, D=distance/abs(DIR(i))
      REAL ROUT(3),RIN(3),DIR(3)
      DO 100 I=1,3
         ROUT(I)=RIN(I)+D*DIR(I)
100   CONTINUE
      RETURN
      END

