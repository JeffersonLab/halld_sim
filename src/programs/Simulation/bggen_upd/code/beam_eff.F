C
      SUBROUTINE BEAM_EFF(IERR)
C
C---   Effective beam spectrum - the beam flux times the cross section
C---   ID - histogram with the dN/dE*sigma(E), 
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
C
      LOGICAL HEXIST
      REAL GPXSECT,GETXSEC,HI,HSUM
      EXTERNAL GPXSECT,GETXSEC,HEXIST,HI,HSUM

      INTEGER i
      REAL flxde,xsec,de,e,xstot,xlum,ecurr,targth,frate,vv
     +    ,fflux
C
      IF(IBEATYP.EQ.0) GO TO 999
C
      IERR=1
C
      IF(.NOT.HEXIST(IDBEAM+11)) THEN
         WRITE(6,1020) IDBEAM+11
 1020    FORMAT('  *** Initialization error: beam flux hist missing',I9)
         GO TO 999
      ELSE
         IF(HSUM(IDBEAM+11).LE.0.) THEN
            WRITE(6,1030) IDBEAM+11,HSUM(IDBEAM+11)
 1030       FORMAT('  *** Initialization error: beam flux hist sum='
     +            ,F10.2)
            GO TO 999
         ENDIF
      ENDIF
C
      CALL HBOOK1(IDBEAM   ,'Beam flux dn/dE*sigma(E)'
     +                  ,NHBEA,EHLIM(1),EHLIM(2),0.)
      CALL HBOOK1(IDBEAM+10,'Cross section in mb'
     +                  ,NHBEA,EHLIM(1),EHLIM(2),0.)

      targth=30.                ! target thickness
      ecurr=350.E-9             ! beam current on a 50um crystal - the real "high luminosity" setting
      xlum=ecurr/1.6E-19*targth*0.071*0.602 ! auxil. luminosity factor:  beam electrons/s * target thickness in 1/b 
C
      de=(EHLIM(2)-EHLIM(1))/NHBEA
C
      fflux=0.
      frate=0.
      DO i=1,NHBEA
         e=EHLIM(1)+(i-0.5)*de
         IF(ISIMUL.EQ.0) THEN
            xstot=GPXSECT(e)*1.E-3   ! x-sec in b  BG
         ELSE
            xstot=GETXSEC(e)         ! x-sec in b - exclusive reaction
         ENDIF
         CALL HF1(IDBEAM+10,e,xstot)
         flxde=HI(IDBEAM+11,i)
         vv=flxde*xstot*xlum
         CALL HF1(IDBEAM,e,vv)
         fflux=fflux+flxde*de
         frate=frate+vv*de
      ENDDO
C
      RATESEC=frate
      WRITE(6,FMT='(//10X,''Rates:'')')
      WRITE(6,1083) ecurr*1.E9,EHLIM,fflux*ecurr/1.6e-19
 1083 FORMAT(10X,'Beam: ',F8.2,' nA e-, gamma in ',2F6.2,' GeV '
     +         ,'  flux ',E11.3,' Hz') 
      WRITE(6,1085) targth
 1085 FORMAT(10X,'Target: ',F8.2,' cm LH2')
      WRITE(6,1090) frate
 1090 FORMAT(10X,'Interaction rate: ',E12.3,' Hz')
C
      CALL HCOPY(IDBEAM,IDBEAM+1,'    ')  ! a copy of the final histogram to be used for HRNDM1
C
      IERR=0
 999  RETURN
      END
