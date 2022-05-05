      SUBROUTINE BG_EVE(IPRI,IERR)
C
C---  Simulates one BG event
C     IPRI>0 - print this event
C
      IMPLICIT NONE
      INTEGER IPRI,IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_proc.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
      INCLUDE 'bg_reac.inc'
      INCLUDE 'bg_statc.inc'
C
      REAL HRNDM1,HXI
      LOGICAL HEXIST
      EXTERNAL HRNDM1,HXI,HEXIST
C
      INTEGER i,j,ip,ier,lout,idnt,lenc(2),nkfsel,ntry
      REAL ebeam,de,ecent
      CHARACTER cent(2)*6,centt(2)*12
C
      REAL ptmp1(4,2),ptmp2(4,MXTRA)   ! auxill. arrays to simplify the HDDM mapping
      INTEGER ifl1(6,2),ifl2(6,MXTRA)
C
C     ------------------------------------------------------------------
C
      lout=6
      IERR=1
      IEVPROC=-1
      INDUM(1)=0
      INDUM(2)=0
      cent(1)='beam  '
      cent(2)='target'
      CPROCNAM='?'
C
      NEVSTAT(1)=NEVSTAT(1)+1
C
C---    Beam energy
C
      IF(IDBEAM.EQ.0.OR..NOT.HEXIST(IDBEAM+1)) GO TO 999
C
      IF(IBEATYP.EQ.0) THEN
         ebeam=EPH_LIM(1)             !   fixed energy
      ELSE
         ntry=0
 100     ntry=ntry+1
         ebeam=HRNDM1(IDBEAM+1)
         IF(ebeam.LT.EPH_LIM(1).OR.ebeam.GT.EPH_LIM(2)) THEN  ! check the boundaries
            IF(ntry.LT.1000) GO TO 100
            WRITE(lout,1000) ebeam,EPH_LIM
 1000       FORMAT(' *** ERROR: Beam energy out of range ',3F9.3)
            GO TO 999
         ENDIF
      ENDIF
C
C---   Beam/target definitions
C
      ITPIN(1,1)=IPINIT(1)    ! beam GEANT type
      ITPIN(1,2)=IPINIT(2)    ! target GEANT type
      ITPIN(2,1)=KF_PART(ITPIN(1,1)) ! KF types
      ITPIN(2,2)=KF_PART(ITPIN(1,2))
C
      DO i=1,2
         AMIN(i)=AM_PART(ITPIN(1,i))
      ENDDO
      DO i=1,3
         PIN(i,1)=0.
         PIN(i,2)=0.
      ENDDO
      PIN(3,1)=ebeam
C
      NTRA=0
C
C---   Choose the package
C
      ierr=0
      IF(ISIMUL.EQ.0) THEN
         IF(ebeam.LT.EPYMIN) THEN
C
            CALL LOWEN_EVE(ierr)
            IF(IEVPROC.GT.0) CPROCNAM=CNPROC(IEVPROC)
C
         ELSE
C
            CALL PYTH_EVE(ierr)
            CPROCNAM='PYTHIA'
C     
         ENDIF
C
      ELSEIF(ISIMUL.EQ.1) THEN
         DO i=1,2
            centt(i)=CH_PART(IPREAC(i))(1:12)
            lenc(i)=LEN_TRIM(centt(i))
         ENDDO
         CPROCNAM='Exclusive production '//centt(1)(1:lenc(1))
     +            //'+'//centt(2)(1:lenc(2))
         CALL REAC_EVE(ierr)
      ENDIF
      IF(ierr.NE.0) GO TO 999
C
      NEVSTAT(2)=NEVSTAT(2)+1
C
C---    Event selection (if needed)
C
      ierr=0
      DO i=1,MXEVSEL
         nkfsel=0
         IF(KFEVSEL(i).NE.0) THEN
            DO ip=1,NTRA
               IF(ITPTRA(3,ip).EQ.KFEVSEL(i)) nkfsel=nkfsel+1 
            ENDDO
            IF(nkfsel.EQ.0) THEN
               ierr=1                ! required particle not found in the record
               GO TO 999
            ENDIF
         ENDIF
      ENDDO
C
      NEVSTAT(3)=NEVSTAT(3)+1
C
C---    Rad corrections
C
      IF(NRADCOR.GT.0) CALL BG_RADCOR
C
C---  Remove the GEANT type for the decaying particles (KF type is retained)
C             needed to avoid copying these particles into GEANT
C
      DO i=1,NTRA
         IF(ITPTRA(2,i).NE.1.AND.ITPTRA(3,i).NE.0) ITPTRA(1,i)=0
      ENDDO  
C
C---     Print the event
C
      IF(IPRI.NE.0) THEN
        CALL BG_PRI(IPRI)
      ENDIF
C
C---   Output file for HDDM
C
      IF(IWROUT(1).NE.0) THEN
         DO i=1,2
            DO j=1,6
               ifl1(j,i)=0
            ENDDO
            ifl1(1,i)=ITPIN(1,i)
            ifl1(3,i)=ITPIN(2,i)
            DO j=1,3
               ptmp1(j,i)=PIN(j,i)
            ENDDO
            ptmp1(4,i)=SQRT(ptmp1(1,i)**2+ptmp1(2,i)**2+ptmp1(3,i)**2
     +                     +AMIN(i)**2)
         ENDDO
         DO i=1,NTRA
            DO j=1,6
               ifl2(j,i)=ITPTRA(j,i)
            ENDDO
C
            DO j=1,3
               ptmp2(j,i)=PTRA(j,i)
            ENDDO
            ptmp2(4,i)=SQRT(ptmp2(1,i)**2+ptmp2(2,i)**2+ptmp2(3,i)**2
     +                     +AMTRA(i)**2)
         ENDDO
         CALL WRITE_HDDM_EVENT(RUNNO, IEVENT,IEVPROC
     +                        ,ifl1(1,1),ptmp1(1,1)
     +                        ,NTRA,ifl2(1,1),ptmp2(1,1))
C         write(6,1010) (i,(ifl1(j,i),j=1,6),(ptmp1(j,i),j=1,4),i=1,2)
C         write(6,1010) (i,(ifl2(j,i),j=1,6),(ptmp2(j,i),j=1,4),i=1,NTRA)
         NEVSTAT(4)=NEVSTAT(4)+1
      ENDIF
C
C---   Sequential output file
C
      IF(IWROUT(2).NE.0) THEN
         WRITE(LUNWR(2)) IEVENT,IEVPROC
     +        ,(( ITPIN(j,i),j=1,2), AMIN(i),( PIN(j,i),j=1,3),i=1,2)
     +   ,NTRA,((ITPTRA(j,i),j=1,6),AMTRA(i),(PTRA(j,i),j=1,3),i=1,NTRA)
      ENDIF
C
C---   NTUPLE
C
      IF(IWROUT(3).NE.0) THEN
         idnt=9
         CALL HFNT(idnt)
         NEVSTAT(5)=NEVSTAT(5)+1
      ENDIF
C
      IERR=0
 999  CONTINUE
C      write(6,*) ebeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C

