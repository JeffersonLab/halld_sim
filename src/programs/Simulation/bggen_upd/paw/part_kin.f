      REAL FUNCTION PART_KIN(IFL,KGEANT,KPYTH,IDH)
C
C--       Fills IDH with the kin. parameters of all tracks of a given type 
C       IFL=0 - p
C          =1 - theta (degrees)
C          =2 - p(Y)-theta(x)(degrees) 
C          =3 - theta rad 
C          =4 - cos(th)  ~ dN/dOmega
C          =4 - th  weight 1/sin(theta)
C       KTYP>0 - GEANT particle type
C          <=0  use KPYTH - PYTHIA KF type
C
      IMPLICIT NONE
      INTEGER IFL,KGEANT,KPYTH,IDH
C
      INCLUDE ?
      LOGICAL HEXIST
C
      INTEGER ip,j,ifirst,ievstart,nfind,ifind
      REAL pf,th,qq,ct,thd
      DATA ifirst/1/
      DATA ievstart/0/
C
      IF(ifirst.EQ.1.OR.IDNEVT.EQ.ievstart) THEN
         IF(IDH.NE.0.AND.HEXIST(IDH)) THEN
            CALL HRESET(IDH,'    ')
            ievstart=IDNEVT
         ELSE
            WRITE(6,*) ' *** ERROR: no histogram ID=',IDH
         ENDIF
      ENDIF
      ifirst=0
C
      nfind=0
C
      DO ip=1,NP
         ifind=0
         IF(KGEANT.GT.0) THEN            
            IF(KGEANT.EQ.ITYP(1,ip)) THEN
               ifind=1
            ENDIF
         ELSE IF(KPYTH.NE.0) THEN
            IF(KPYTH.EQ.ITYP(3,ip)) THEN
               ifind=1
            ENDIF
         ENDIF
C         write(6,*) ifind,KGEANT,KPYTH,ITYP(1,ip),ITYP(3,ip)
         IF(ifind.NE.0) THEN
            nfind=nfind+1
            qq=0.
            DO j=1,3
               qq=qq+POUT(j,ip)**2
            ENDDO
            pf=SQRT(qq)
            ct=POUT(3,ip)/pf
            th=ACOS(ct)
            thd=th*180./3.1416
C
            IF(IFL.EQ.0) THEN
               CALL HFILL(IDH,pf,0.,1.)
            ELSE IF(IFL.EQ.1) THEN
               CALL HFILL(IDH,thd,0.,1.)
            ELSE IF(IFL.EQ.2) THEN
               CALL HFILL(IDH,th,pf,1.)
            ELSE IF(IFL.EQ.3) THEN
               CALL HFILL(IDH,th,0.,1.)
            ELSE IF(IFL.EQ.4) THEN
               CALL HFILL(IDH,ct,0.,1.)
            ELSE IF(IFL.EQ.5) THEN
               IF(th.GT.0.) CALL HFILL(IDH,th,0.,1./SIN(th))
            ENDIF
         ENDIF
      ENDDO
C
      PART_KIN=nfind
C
      RETURN
      END
