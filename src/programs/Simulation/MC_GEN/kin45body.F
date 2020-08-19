c******************************************************************************
c     KIN45BODY.F
c
c     This is a wrapper for the CERNLIB program GENBOD.F 
c     It computes 4-body or 5-body phase space using GENBOD for the kinematic 
c     calculations.
c
c     GENBOD is included below, with partial small modifications to modernize it.
c     It can be linked in to the code or not, depending on whether you want
c     to use CERNLIB directly.
c     
c     R. A. Schumacher, CMU 5-08
c
c     5-15-15 - make the 4-vectors have 5 elements, with E in element 5
c     2-11-18 - fix a bug in the way the total energy was being transferred to local variables
c     8-01-20 - use the CERNLIB GENBOD code that is included below.
c
      subroutine kin45body(nbody,p0,p1,p2,p3,p4,p5,ierr)
      real p0(5),p1(5),p2(5),p3(5),p4(5),p5(5),beta(3)
      real p1tmp(5),p2tmp(5),p3tmp(5),p4tmp(5),p5tmp(5)
      data valmax/0/
      common/kin4bodylocal/valmax
      common/monte/iseed
      COMMON/GENIN / NP, TECM, AMASS(18), KGENEV
      COMMON/GENOUT/ PCM(5,18) , WT
c     
      x = p0(1)   ! dummy to fool compiler
      kgenev = 1  ! "constant cross section with energy"
      np = nbody  ! define number of particles to generate
      amass(1) = p1(4)  !load the masses for genbod call
      amass(2) = p2(4)
      amass(3) = p3(4)
      amass(4) = p4(4)
      amass(5) = p5(4)
      p1tmp(4) = p1(4)
      p2tmp(4) = p2(4)
      p3tmp(4) = p3(4)
      p4tmp(4) = p4(4)
      p5tmp(4) = p5(4)
      tecm     = p0(4) !decaying state at rest has the mass of p0
c
c     Double check of mass limits
c
      rmtot = p1(4)+p2(4)+p3(4)+p4(4)+p5(4)
      if(rmtot.ge.p0(4))then
c         write(6,*)'Bloody murder',rmtot,p0(4),p1(4),p2(4),p3(4),p4(4),p5(4)
         ierr = 1
         return
      endif
c
      call genbod        !computes phase space decay in rest frame of p0
c      write(6,*)'Phase space GENBOD calculation not linked into program'
c      call exit()
c
      do i=1,3
         p1tmp(i) = pcm(i,1)
         p2tmp(i) = pcm(i,2)
         p3tmp(i) = pcm(i,3)
         p4tmp(i) = pcm(i,4)
         p5tmp(i) = pcm(i,5)
      enddo
      p1tmp(5) = pcm(4,1)  !our total energy E is GENBOD's item 4
      p2tmp(5) = pcm(4,2)
      p3tmp(5) = pcm(4,3)
      p4tmp(5) = pcm(4,4)
      p5tmp(5) = pcm(4,5)
c
c     Boost all particles to the frame wherein the original particle, p0, 
c     is moving.  Start by computing the beta of the decaying particle
c     that defined the boost velocity
c
      energy = etot(p0)  ! total energy of decaying particle
      do i=1,3
         beta(i) = -p0(i)/energy
      enddo
      call boost(p1tmp,beta,p1)
      call boost(p2tmp,beta,p2)
      call boost(p3tmp,beta,p3)
      call boost(p4tmp,beta,p4)
      if(nbody.eq.5)call boost(p5tmp,beta,p5)
c
c     Kinematic checks...momentum conservation, etc.
c     
c      write(6,*)'p0:',p0, wt
c      write(6,*)'p1:',p1,p1tmp
c      write(6,*)'p2:',p2,p2tmp
c      write(6,*)'p3:',p3,p3tmp
c      write(6,*)'p4:',p4,p4tmp
c      write(6,*)'p5:',p5,p5tmp
c      write(6,*)'Beta',beta
c      dpx = p1(1)+p2(1)+p3(1)+p4(1)+p5(1) - p0(1)
c      dpy = p1(2)+p2(2)+p3(2)+p4(2)+p5(2) - p0(2)
c      dpz = p1(3)+p2(3)+p3(3)+p4(3)+p5(3) - p0(3)
c      dpxtmp = p1tmp(1)+p2tmp(1)+p3tmp(1)+p4tmp(1)+p5tmp(1) 
c      dpytmp = p1tmp(2)+p2tmp(2)+p3tmp(2)+p4tmp(2)+p5tmp(2) 
c      dpztmp = p1tmp(3)+p2tmp(3)+p3tmp(3)+p4tmp(3)+p5tmp(3) 
c      write(6,*)'Momentum check rest   frame :',dpxtmp,dpytmp,dpztmp
c      write(6,*)'Momentum check moving frame :',dpx,dpy,dpz
c      write(6,20)'Genbod weight',wt
c 20   format(1x,a,f8.4)
      ierr = 0
      return
      end 

c******************************************************************************      
*
* $Id: genbod.F,v 1.1.1.1 1996/03/22 16:42:49 mclareni Exp $
*
* $Log: genbod.F,v $
* Revision 1.1.1.1  1996/03/22 16:42:49  mclareni
* Phtools
*
*
      SUBROUTINE GENBOD         ! use this code
c     SUBROUTINE xGENBOD   ! use CERNLIB's version
C   SUBROUTINE TO GENERATE N-BODY EVENT
C   ACCORDING TO FERMI LORENTZ-INVARIANT PHASE SPACE
C   ADAPTED FROM FOWL (CERN W505) SEPT. 1974 BY F. JAMES
C   EVENTS ARE GENERATED IN THEIR OWN CENTER-OF-MASS,
C   BUT MAY BE TRANSFORMED TO ANY FRAME USING LOREN4
C
C   INPUT TO SUBROUTINE IS THRU COMMON BLOCK GENIN
C             NP=NUMBER OF OUTGOING PARTICLES (.LT. 19)
C             TECM=TOTAL ENERGY IN CENTER-OF-MASS
C             AMASS(I)=MASS OF ITH OUTGOING PARTICLE
C             KGENEV=1 FOR CONSTANT CROSS SECTION
C                      2 FOR FERMI ENERGY-DEPENDANCE
C
C   OUTPUT FROM SUBROUTINE IS THRU COMMON BLOCK GENOUT
C             PCM(1,I)=X-MOMENTUM IF ITH PARTICLE
C             PCM(2,I)=Y-MOMENTUM IF ITH PARTICLE
C             PCM(3,I)=Z-MOMENTUM IF ITH PARTICLE
C             PCM(4,I)=ENERGY OF ITH PARTICLE
C             PCM(5,I)=MOMENTUM OF ITH PARTICLE
C             WT=WEIGHT OF EVENT
      COMMON/GENIN / NP, TECM, AMASS(18), KGENEV
      COMMON/GENOUT/ PCM(5,18) , WT
      DIMENSION EMM(18)
      DIMENSION RNO(50)
C--       PCM1 IS LINEAR EQUIV. OF PCM TO AVOID DOUBLE INDICES
      DIMENSION EM(18),PD(18),EMS(18),SM(18),FFQ(18),PCM1(90)
      EQUIVALENCE (NT,NP),(AMASS(1),EM(1)),(PCM1(1),PCM(1,1))
C FFQ(N)=PI * (TWOPI)**(N-2) / (N-2)FACTORIAL
      DATA FFQ/0.,3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
     2                       256.3704, 268.4705, 240.9780, 189.2637,
     3                       132.1308,  83.0202,  47.4210,  24.8295,
     4                        12.0006,   5.3858,   2.2560,   0.8859/
      DATA KNT,TWOPI/0,6.2831853073/
C        INITIALIZATION
      KNT=KNT + 1
      IF(KNT.GT.1) GOTO 100
      WRITE(6,1160)
      WRITE(6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
 100  CONTINUE
      IF(NT.LT.2) GOTO 1001
      IF(NT.GT.18) GOTO 1002
      NTM1=NT-1
      NTM2=NT-2
      NTP1=NT+1
      NTNM4=3*NT - 4
      EMM(1)=EM(1)
      TM=0.0
      DO 200 I=1,NT
         EMS(I)=EM(I)**2
         TM=TM+EM(I)
         SM(I)=TM
 200  continue
C     CONSTANTS DEPENDING ON TECM
      TECMTM=TECM-TM
      IF(TECMTM.LE.0.0) GOTO 1000
      EMM(NT)=TECM
      IF(KGENEV.eq.1) then
C        CONSTANT CROSS SECTION AS FUNCTION OF TECM
         EMMAX=TECMTM+EM(1)
         EMMIN=0.0
         WTMAX=1.0
         DO 350 I=2,NT
            EMMIN=EMMIN+EM(I-1)
            EMMAX=EMMAX+EM(I)
            WTMAX=WTMAX* xPDK(EMMAX,EMMIN,EM(I))
 350     continue
         WTMAXQ=1.0/WTMAX
      else
C--      FERMI ENERGY DEPENDENCE FOR CROSS SECTION
         WTMAXQ=TECMTM**NTM2*FFQ(NT) / TECM
C        CALCULATION OF WT BASED ON EFFECTIVE MASSES EMM
      endif
c
C--   FILL RNO WITH 3*NT-4 RANDOM NUMBERS,
C--   OF WHICH THE FIRST NT-2 ARE ORDERED.
      iseed = 0
      DO I= 1, NTNM4
c         RNO(I)=RNDM(I)
         RNO(I)=ran(iseed) 
      enddo
      if(ntm2.lt.0)return
      if(ntm2.gt.0)then
         CALL FLPSOR(RNO,NTM2)
         DO J=2,NTM1
            EMM(J)=RNO(J-1)*(TECMTM)+SM(J)
         enddo
      endif
      WT=WTMAXQ
      IR=NTM2
      DO 530 I=1,NTM1
      PD(I)= xPDK(EMM(I+1),EMM(I),EM(I+1))
  530 WT=WT*PD(I)
C--       COMPLETE SPECIFICATION OF EVENT (RAUBOLD-LYNCH METHOD)
      PCM(1,1)=0.0
      PCM(2,1)=PD(1)
      PCM(3,1)=0.0
      DO 570 I=2,NT
      PCM(1,I)=0.0
      PCM(2,I)=-PD(I-1)
      PCM(3,I)=0.0
      IR=IR+1
      BANG=TWOPI*RNO(IR)
      CB=COS(BANG)
      SB=SIN(BANG)
      IR=IR+1
      C=2.0*RNO(IR)-1.0
      S=SQRT(1.0-C*C)
      IF(I.EQ.NT) GOTO 1567
      ESYS=SQRT(PD(I)**2+EMM(I)**2)
      BETA=PD(I)/ESYS
      GAMA=ESYS/EMM(I)
      DO 568 J=1,I
      NDX=5*J - 5
      AA= PCM1(NDX+1)**2 + PCM1(NDX+2)**2 + PCM1(NDX+3)**2
      PCM1(NDX+5)=SQRT(AA)
      PCM1(NDX+4)=SQRT(AA+EMS(J))
      CALL xROTES2(C,S,CB,SB,PCM,J)
      PSAVE=GAMA*(PCM(2,J)+BETA*PCM(4,J))
  568 PCM(2,J)=PSAVE
      GOTO 570
 1567 DO 1568 J=1,I
      AA=PCM(1,J)**2 + PCM(2,J)**2 + PCM(3,J)**2
      PCM(5,J)=SQRT(AA)
      PCM(4,J)=SQRT(AA+EMS(J))
      CALL xROTES2(C,S,CB,SB,PCM,J)
 1568 CONTINUE
  570 CONTINUE
      RETURN
C          ERROR RETURNS
 1000 WRITE(6,1100)
      GOTO 1050
 1001 WRITE(6,1101)
      GOTO 1050
 1002 WRITE(6,1102)
 1050 WRITE(6,1150) KNT
      WRITE(6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
      STOP
 1100 FORMAT(28H0 AVAILABLE ENERGY NEGATIVE )
 1101 FORMAT(33H0 LESS THAN 2 OUTGOING PARTICLES )
 1102 FORMAT(34H0 MORE THAN 18 OUTGOING PARTICLES )
 1150 FORMAT(47H0 ABOVE ERROR DETECTED IN GENBOD AT CALL NUMBER,I7)
 1160 FORMAT(34H0 FIRST CALL TO SUBROUTINE GENBOD )
 1200 FORMAT(36H  INPUT DATA TO GENBOD.         NP=   ,I6/
     +  ,8H   TECM=,E16.7,18H  PARTICLE MASSES=,5E15.5/(42X,5E15.5)
     +)
      END

*
* $Id: pdk.F,v 1.1.1.1 1996/03/22 16:42:47 mclareni Exp $
*
* $Log: pdk.F,v $
* Revision 1.1.1.1  1996/03/22 16:42:47  mclareni
* Phtools
*
*
c      FUNCTION PDK(A,B,C)
      FUNCTION xPDK(A,B,C)
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C--  CALLED FROM -  GENEV
C     PDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
      A2 = A*A
      B2 = B*B
      C2 = C*C
      xPDK = 0.5*SQRT(A2 + (B2-C2)**2/A2 - 2.0*(B2+C2))
      RETURN
      END


*
* $Id: rotes2.F,v 1.1.1.1 1996/03/22 16:42:47 mclareni Exp $
*
* $Log: rotes2.F,v $
* Revision 1.1.1.1  1996/03/22 16:42:47  mclareni
* Phtools
*
*
c      SUBROUTINE ROTES2(C,S,C2,S2,PR,I)
      SUBROUTINE xROTES2(C,S,C2,S2,PR,I)
C--  CALLED FROM - GENEV
C         THIS SUBROUTINE NOW DOES TWO ROTATIONS (XY AND XZ)
      DIMENSION PR(50)
      K1 = 5*I - 4
      K2 = K1 + 1
      SA = PR(K1)
      SB = PR(K2)
      A      = SA*C - SB*S
      PR(K2) = SA*S + SB*C
      K2 = K2 + 1
      B = PR(K2)
      PR(K1) = A*C2 - B*S2
      PR(K2) = A*S2 + B*C2
      RETURN
      END


* $Log: flpsor.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
c#include "kerngen/pilot.h"
      SUBROUTINE FLPSOR(A,N)
C
C CERN PROGLIB# M103    FLPSOR          .VERSION KERNFOR  3.15  820113
C ORIG. 29/04/78
C
C   SORT THE ONE-DIMENSIONAL FLOATING POINT ARRAY A(1),...,A(N) BY
C   INCREASING VALUES
C
C-    PROGRAM  M103  TAKEN FROM CERN PROGRAM LIBRARY,  29-APR-78
C
      DIMENSION A(N)
      COMMON /SLATE/ LT(20),RT(20)
      INTEGER R,RT
C
      LEVEL=1
      LT(1)=1
      RT(1)=N
   10 L=LT(LEVEL)
      R=RT(LEVEL)
      LEVEL=LEVEL-1
   20 IF(R.GT.L) GO TO 200
      IF(LEVEL) 50,50,10
C
C   SUBDIVIDE THE INTERVAL L,R
C     L : LOWER LIMIT OF THE INTERVAL (INPUT)
C     R : UPPER LIMIT OF THE INTERVAL (INPUT)
C     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)
C     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)
C
  200 I=L
      J=R
      M=(L+R)/2
      X=A(M)
  220 IF(A(I).GE.X) GO TO 230
      I=I+1
      GO TO 220
  230 IF(A(J).LE.X) GO TO 231
      J=J-1
      GO TO 230
C
  231 IF(I.GT.J) GO TO 232
      W=A(I)
      A(I)=A(J)
      A(J)=W
      I=I+1
      J=J-1
      IF(I.LE.J) GO TO 220
C
  232 LEVEL=LEVEL+1
      IF((R-I).GE.(J-L)) GO TO 30
      LT(LEVEL)=L
      RT(LEVEL)=J
      L=I
      GO TO 20
   30 LT(LEVEL)=I
      RT(LEVEL)=R
      R=J
      GO TO 20
   50 RETURN
      END

c******************************************************************************
