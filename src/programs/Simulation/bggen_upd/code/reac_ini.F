      SUBROUTINE REAC_INI(IERR)
C
C---   Exclusive reaction initialization
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_reac.inc'
C
      REAL GPTDIST
C
      INTEGER i,j,nb,idt,lout,ier,ir,ntrng,itsng
      REAL e0,e1,e0c,t1,t2,tt,dt,pdfu
C
      CHARACTER crng*2,clim*6
C
      IERR=1
      lout=6
C
      ier=0
      DO i=1,2
         IF(IPREAC(i).LT.1.OR.IPREAC(i).GT.MXPART) THEN
            ier=1
         ELSE
            IF(KF_PART(IPREAC(i)).EQ.0) ier=1
         ENDIF
      ENDDO
      IF(ELREAC(1).GT.EPH_LIM(2).OR.ELREAC(2).LT.EPH_LIM(1).OR.
     +   ELREAC(1).LT.0..OR.ELREAC(2).LT.0.) ier=2 
      IF(NPXREAC.LT.2) ier=3
      IF(TLREAC(1).GE.TLREAC(2)) ier=4
      ntrng=0                   ! number of active t ranges
      itsng=0                   ! set for a full  active range (covering all energies) if only the first 2 coefficients are set (a simple exponent)
      e0=MAX(ELREAC(1),EPH_LIM(1)) ! lower beam energy limit
      e1=MIN(ELREAC(2),EPH_LIM(2)) ! upper beam energy limit
      e0c=e0                    ! the current upper limit
      DO ir=1,MXTRNG
         IF(TPAREAC(1,ir).GT.e0c) THEN ! inside the useful energy range
            IF(TPAREAC(2,ir).GT.0.) THEN ! normaliation constant non zero
               ITDRNG(ir)=1
               IF(itsng.EQ.0) THEN
                  itsng=ir
                  DO j=4,11
                     IF(ABS(TPAREAC(j,ir)).GT.1.E-8) itsng=0
                  ENDDO
                  IF(itsng.GT.0) THEN
                     ITDRNG(ir)=2
                  ENDIF
               ENDIF
               ntrng=ntrng+1
               e0c=TPAREAC(1,ir)
            ENDIF
         ENDIF
         if(e0c.GE.e1) GO TO 30
      ENDDO
 30   IF(ntrng.EQ.0) ier=5
      WRITE(lout,1025) TLREAC
 1025 FORMAT(' t-distributions. Exclusive reaction: t limits ',2F10.3
     +        ,/,'  #  use Emax_beam   parameters ')   
      DO ir=1,MXTRNG
         IF(ITDRNG(ir).EQ.2.AND.ntrng.GT.1) ITDRNG(ir)=1
         WRITE(lout,1026) ir,ITDRNG(ir),(TPAREAC(j,ir),j=1,11)
 1026    FORMAT(I3,1X,I2,2X,F7.3,3X,10F10.3)
      ENDDO
      IF(ier.NE.0) THEN
         WRITE(6,1030) IPREAC,TLREAC,ELREAC,NPXREAC
 1030    FORMAT('  *** Initialization error REAC: :',
     +                                  4I4,3E11.3,I6)
         GO TO 999
      ENDIF
      DO i=1,NPXREAC
         ESREAC(i)=ELREAC(1)+(ELREAC(2)-ELREAC(1))/(NPXREAC-1)*(i-1)
C         write(6,*) i,ESREAC(i)
      ENDDO
C
      ITDSING=0
      TSLREAC=0.
      DO ir=1,MXTRNG
C---            A simple single exponential function, no histogram is needed
         IF(ITDRNG(ir).EQ.2) THEN
            ITDSING=ir
            TSLREAC=TPAREAC(2,ir)
C
C---            A complex function for the t distribution, use HBOOK 
         ELSEIF(ITDRNG(ir).EQ.1) THEN
            idt=IDTDIST+10*ir
            nb=1000
            t1=TLREAC(1)
            t2=TLREAC(2)
            dt=(t2-t1)/nb
            WRITE(crng,FMT='(I2)') ir
            WRITE(clim,FMT='(F6.3)') TPAREAC(1,ir)
            CALL HBOOK1(idt
     +       ,'t-distribution, range '//crng//' Eb<'//clim,nb,t1,t2,0.)
            tt=t1
            DO i=1,nb
               tt=t1+(i-0.5)*dt
               pdfu=GPTDIST(tt,TPAREAC(2,ir),10)   ! t-distribution function
               CALL HF1(idt,tt,pdfu)
            ENDDO
         ENDIF
      ENDDO
C
      IERR=0
 999  RETURN
      END
