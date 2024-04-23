C
      SUBROUTINE BEAM_INI(IERR)
C
C---   Sets the effective beam energy spectrum, as a product of the
C      beam energy spectrum and the cross section of the reactions to study 
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_beamext.inc'
C
      REAL DNIDX,DNCDX,DIVDIF
      EXTERNAL DNIDX,DNCDX,DIVDIF

      INTEGER i,ipol
      REAL flxdp,xsec,dx,p,dp
C
      IERR=1
C
      IF(EPH_LIM(1).GT.EPH_LIM(2)) THEN
         WRITE(6,1005) EPH_LIM
 1005    FORMAT('  *** Initialization error: energy limits:',2F10.4)
         GO TO 999
      ELSE IF(EPH_LIM(1).EQ.EPH_LIM(2)) THEN
         IBEATYP=0
C
C---      Increase E2 slightly in order to make a valid histogram
C
         EPH_LIM(2)=EPH_LIM(1)*1.0001
C
      ELSE
         IF(NBMEXT.LE.0) THEN          !  Bresstrahlung beam
            IBEATYP=1
C
C---       Bremsstrahlung beam: the E0 and Epeak should be consistent
C
            IF(EELEC.LT.EPH_LIM(2)) THEN
               WRITE(6,1010) EELEC,EPH_LIM(2)
 1010          FORMAT('  *** Initialization error: Ee<E2:',2F10.4)
               GO TO 999
            ENDIF
C
            IF(EPEAK.GE.EELEC) THEN
               WRITE(6,1020) EPEAK,EELEC
 1020          FORMAT('  *** Warning: peak>Ee:',2F10.4
     +              ,' no coherent Bremmstrahlung')
            ENDIF
C
         ELSE       ! user-defined beam spectrum
            IF(NBMEXT.LT.2) THEN
               WRITE(6,1040) NBMEXT
 1040          FORMAT('  *** Initialization error: NBEAMEXT='
     +              ,I5,'<2')
               GO TO 999
            ENDIF
            IBEATYP=2
            IF(IBMEXT.LT.0.OR.IBMEXT.GT.2) THEN
               WRITE(6,1048) IBMEXT
 1048          FORMAT('  *** Initialization error: flag '
     +               ,'IBMEXT=',I4,' unknown')
               GO TO 999
            ELSE IF(IBMEXT.GT.0) THEN
               EXT_LIM(1)=PBEAEXT(1)
               EXT_LIM(2)=PBEAEXT(NBMEXT)
               WRITE(6,*) ' --- Beam range EXT_LIM are set to ',EXT_LIM 
            ENDIF
            IF(EXT_LIM(1).GE.EXT_LIM(2)) THEN
               WRITE(6,1045) EXT_LIM
 1045          FORMAT('  *** Initialization error: energy limits '
     +               ,'BEAMEXT ',2F10.4)
               GO TO 999
            ENDIF
            IF(EXT_LIM(1).GT.EPH_LIM(1).OR.
     +         EXT_LIM(2).LT.EPH_LIM(2)) THEN
               WRITE(6,1046) EXT_LIM,EPH_LIM
 1046          FORMAT('  *** Initialization error: energy limits '
     +               ,'BEAMEXT should cover the limits EPHLIM',4F10.4)
               GO TO 999
            ENDIF
         ENDIF
C
      ENDIF
C
C---         For the background simulation: check the lower limit
C
      IF(ISIMUL.EQ.0) THEN
         IF(EPH_LIM(1).LT.ELOWMIN) THEN
            WRITE(6,1050) EPH_LIM(1),ELOWMIN
 1050       FORMAT('  --- Initialization error: EPH_LIM(1)'
     +            ,' < ELOWMIN:',2F10.4)
            GO TO 999
         ENDIF
      ENDIF
C
C---    The range and the number of bins
C
      DO i=1,2
         PHLIM(i)=EPH_LIM(i)
      ENDDO
      dp=MIN(0.005,PHLIM(2)-PHLIM(1))
      NHBEA=INT((PHLIM(2)-PHLIM(1))/dp+0.001) ! number of bins in the energy tables
C
C--- Exceptions
C
      IF(IBEATYP.EQ.0) THEN
         NHBEA=1   !   One bin - fixed energy at EPH_LIM(1)
      ELSE IF(IBEATYP.EQ.2) THEN                         ! Beam spectrum set by the user
         IF(IBMEXT.EQ.0) THEN
            DO i=1,2
               PHLIM(i)=EXT_LIM(i)
            ENDDO
            NHBEA=NBMEXT
            dp=(PHLIM(2)-PHLIM(1))/NHBEA
         ENDIF
      ENDIF
C
      IF(IDBEAM.EQ.0) THEN
         WRITE(6,1060) IDBEAM
 1060    FORMAT('  *** Initialization error: IDBEAM =',I3)
         GO TO 999
      ENDIF
         
      IF(IBEATYP.GT.0) THEN
         CALL HBOOK1(IDBEAM+11
     +              ,'Beam flux dn/dp' ,NHBEA,PHLIM(1),PHLIM(2),0.)
C
         IF(IBEATYP.EQ.1) THEN
C---   Initialize the coherent Bremsstrahlung 
            CALL COBREMS(EELEC,EPEAK,EEMIT,RADT,ZCOLL/100.,DCOLL,0) ! collimator distance in m
         ENDIF
C
         dx=dp/EELEC
         DO i=1,NHBEA
            p=PHLIM(1)+(i-0.5)*dp
            IF(IBEATYP.EQ.1) THEN
               flxdp=(DNIDX(p/EELEC)+DNCDX(p/EELEC))*dx/dp ! dN/dx - photon rate per one incident electron (E=p)
            ELSE IF(IBEATYP.EQ.2) THEN
               IF(IBMEXT.EQ.0) THEN
                  flxdp=RBEAEXT(i)     ! The bin contents
               ELSE
                  ipol=3
                  flxdp=DIVDIF(RBEAEXT(1),PBEAEXT(1),NBMEXT,p,ipol) ! interpolation between points
               ENDIF
            ENDIF
            CALL HF1(IDBEAM+11,p,flxdp)
         ENDDO
C
         CALL HCOPY(IDBEAM+11,IDBEAM+12,'    ')  ! for HRNDM
      ENDIF

      IERR=0
 999  RETURN
      END
