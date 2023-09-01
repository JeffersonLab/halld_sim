C
C---  Simulates "background" photoproduction by a coherent Bremsstrahlung beam
C     Reaction: gamma+p
C     Control flags are read from fort.15 (FFREAD)
C     E>3 GeV (can be changed in fort.15) - use PYTHIA
C     E<3 GeV - a coctail of several dominating photoproduction processes        
C
C     Includes: a) calculation of the coherent+incoherent photon energy spectra
C               b) parametrization for the cross sections:
C                  - full (formula fit to data)
C                  - p pi0, n pi+ - using SAID
C                  - p 2pi, n 2pi, p eta, p 3pi, n 3pi (formula fit to data)
C               c)  simulation if unbiased (equal weight) events in a given beam energy range: 
C                  - beam energy simulated (beam spectrum times the total cross section)
C                  - the process is chosen randomly accordingly to the their probabilities 
C
C---  Input: file "fort.15" 
C                      the number of events, the beam energy range, the distance to the collimator etc
C            file "particle.dat" contains a table for particle masses (GEANT numbering)
C
C
      SUBROUTINE BGGEN
C
      IMPLICIT NONE
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_proc.inc'
      INCLUDE 'bg_evec.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_statc.inc'
C
      INTEGER ierr,iev,ipri,i,kc,kcs,nevsel
      INTEGER nproc(0:MXPROC)
      CHARACTER cnam(0:MXPROC)*16,cmom*18,cpnam*16
      REAL    rat,ama,wgt
      INTEGER kct,ndc,ifd
C
C     ------------------------------------------------------------------
C
      CALL BG_INI(ierr)
      IF(ierr.NE.0) GO TO 999
C
      DO i=0,MXPROC
         nproc(i)=0
      ENDDO
      IEVENT=0
C
      DO iev=1,NEVENT
C
         IEVENT=IEVENT+1
         ipri=0
         IF(iev.LE.NPRIEV) ipri=1
C
         CALL BG_EVE(ipri,ierr)
C
         IF(ierr.EQ.0) THEN
            IF(IEVPROC.GE.0.AND.IEVPROC.LE.MXPROC) THEN
               nproc(IEVPROC)=nproc(IEVPROC)+1
            ENDIF
         ELSE
            IEVENT=IEVENT-1  ! the event is rejected
         ENDIF
C
      ENDDO
C
      WRITE(6,1980) 
 1980 FORMAT(///1X,100('='))
      WRITE(6,1990) (NEVSTAT(i),i=1,5),RATESEC/1000.
 1990 FORMAT(
     +     3X,I12,' events initialy simulated in the full energy range'
     +  ,/,3X,I12,' events successfully simulated'
     +  ,/,3X,I12,' events selected'
     +  ,/,3X,I12,' events written out to HDDM'
     +  ,/,3X,I12,' events written out to HBOOK NTuple'
     +  ,/,3X,F12.2,' reference interaction rate, kHz')
      IF(IEVENT.GT.0) THEN
         cnam(0)='PYTHIA '
         DO i=1,MXPROC
            cnam(i)=CNPROC(i)
         ENDDO
         WRITE(6,2000) 
 2000    FORMAT('  process ',16X,'     events   fraction   range')
         WRITE(6,2005) 
 2005    FORMAT(1X,100('-'))
         DO i=0,MXPROC
            IF(nproc(i).GT.0) THEN
               cmom=' '
               IF(i.EQ.0) THEN
                  WRITE(cmom,2006) MAX(EPYMIN,EPH_LIM(1)),EPH_LIM(2)
 2006             FORMAT(F5.2,'<E<',F5.2,' GeV')
               ELSE IF(i.EQ.1) THEN
                  WRITE(6,2005) 
                  WRITE(cmom,2006) EPH_LIM(1),EPYMIN
               ENDIF
               WRITE(6,2010) i,cnam(i),nproc(i)
     +                ,nproc(i)/REAL(IEVENT)*100.,cmom
            ENDIF
         ENDDO
 2010    FORMAT(3X,I4,2X,A16,2X,I10,3X,F5.1,' %',2X,A)
      ENDIF
      IF(nproc(0).GT.0) THEN
         WRITE(6,2005)
         WRITE(6,2012)
 2012    FORMAT(' PYTHIA total multiplicity in simulated events'
     +    ,/,'      multiplicity          kf   %undecayed   kc'
     +    ,'   mass    name ')     
         DO kc=1,MXPKC
            kcs=kc
            DO i=1,2
               IF(MULT_KC(kcs).GT.0) THEN
                 rat=NODC_KC(kcs)/REAL(MULT_KC(kcs))*100.
                 CALL PYTH_GETP(KF_KC(kcs),kct,cpnam,ama,wgt,ndc,ifd)
                 WRITE(6,2015) MULT_KC(kcs),KF_KC(kcs),rat,kc,ama,cpnam
 2015            FORMAT(6X,I11,3X,I10,3X,F8.3,3X,I4,1X,F8.4,2X,A)
                 kcs=kcs*(-1)
              ENDIF
           ENDDO
        ENDDO
      ENDIF
      WRITE(6,2020) 
 2020 FORMAT(1X,100('-')///)
C
      CALL BG_END
C
 999  CONTINUE
C
      END

