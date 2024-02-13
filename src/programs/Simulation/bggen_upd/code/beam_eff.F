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
      INCLUDE 'bg_partc.inc'
C
      LOGICAL HEXIST
      REAL GPXSECT,GETXSEC,HI,HSUM
      EXTERNAL GPXSECT,GETXSEC,HEXIST,HI,HSUM

      INTEGER i,nw
      REAL flxde,xsec,de,e,xstot,xlum,ecurr,targth,frate,vv
     +    ,fflux
     +    ,w1,w2,dw,ww,pp
     +    ,ami(2)     ! masses of the beam, terget particles
     +    ,xsmax
      CHARACTER tit*80
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
      
C
      tit='Beam flux dn/dp*sigma(p)'
      CALL HBOOK1(IDBEAM   ,tit,NHBEA,PHLIM(1),PHLIM(2),0.)
      tit='Cross section in mb'
      CALL HBOOK1(IDBEAM+10,tit,NHBEA,PHLIM(1),PHLIM(2),0.)

      targth=30.                ! target thickness
      ecurr=350.E-9             ! beam current on a 50um crystal - the real "high luminosity" setting
      xlum=ecurr/1.6E-19*targth*0.071*0.602 ! auxil. luminosity factor:  beam electrons/s * target thickness in 1/b 
C
      de=(PHLIM(2)-PHLIM(1))/NHBEA
C
      fflux=0.
      frate=0.
      DO i=1,NHBEA
         e=PHLIM(1)+(i-0.5)*de
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
C---   xsec(W) used for reactions on "Fermi-moving" nucleons 
C
      IF(NUCLTARG.GT.0) THEN
C---          xsec(W) in a very wide range
         w1=0.
         w2=15. ! arbitrary: equivalent to a 200GeV beam
         dw=0.001 ! 1MeV
         nw=(w2-w1)/dw
         CALL HBOOK1(IDBEAM+13,'Xsec W dependence',nw,w1,w2,0.)
         ami(1)=AM_PART(IPINIT(1))
         ami(2)=AM_PART(IPINIT(2))
         xsmax=0.
         DO i=1,nw
            ww=w1+dw*(i-0.5)
            e=(ww**2-ami(1)**2-ami(2)**2)/2./ami(2)  !  energy of the beam assuming the target is at rest
            pp=SQRT(e**2-ami(1)**2) ! momentum
C
            IF(ISIMUL.EQ.0) THEN
               xstot=GPXSECT(pp)*1.E-3 ! x-sec
            ELSE
               xstot=GETXSEC(pp) ! x-sec - exclusive reaction
            ENDIF
            xsmax=MAX(xsmax,xstot)
C
            CALL HF1(IDBEAM+13,ww,xstot)
         ENDDO
         IF(xsmax.GT.0.) THEN
            CALL HOPERA(IDBEAM+13,'+',IDBEAM+13,IDBEAM+13,1./xsmax,0.)
         ENDIF
      ENDIF
C
      RATESEC=frate
      WRITE(6,FMT='(//10X,''Rates:'')')
      WRITE(6,1083) ecurr*1.E9,PHLIM,fflux*ecurr/1.6e-19
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
