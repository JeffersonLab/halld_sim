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
      REAL HRNDM1,HX,RNDM
      LOGICAL HEXIST
      EXTERNAL HRNDM1,HX,HEXIST,RNDM
C
      INTEGER i,j,ip,lout,idnt,lenc(2),nkfsel,ntry0,ntry,npfinal
     +       ,itransl   ! >0 - for "Fermi motion" the frame changed
      REAL pbeam,thdg,ptot,pt2
      CHARACTER cent(2)*6,centt(2)*12
C
      REAL pilab(4,2)   ! beam and target in Lab (for Fermi motion)
     +    ,bettar(4)    ! beta of the target rest frame
     +    ,pitar(4,2)   ! beam and target in the target rest frame
     +    ,pitarr(4,2)  ! beam and target in the target rest frame beam along Z
     +    ,rotz(3,3)    ! rotation matrix tar -> tar, beam along Z
     +    ,psf,esf      ! generated (spectral func): nucleon momentum and the removal energy
     +    ,eneff,ameff  ! energy and effective mass of the struck nucleon
     +    ,ct,st,phi,twopi       ! angles
     +    ,xsec         ! cross section
     +    ,ww           ! CM energy
     +    ,xsrnd        ! random value
     +    ,potarr(4)    ! auxil 
     +    ,potar(4)     ! auxil 
     +    ,polab(4)     ! auxil 
     +    ,amtarini     ! the mass of the target nucleon on-shell 
C
      REAL ptmp1(4,2),ptmp2(4,MXTRA)   ! auxill. arrays to simplify the HDDM mapping
      INTEGER ifl1(6,2),ifl2(6,MXTRA)
C
      CHARACTER tit*80
      REAL p1,p2,e1,e2          ! Parameters of the IDNUCLT fine histogram
      SAVE p1,p2,e1,e2
      INTEGER np,ne,nwt,lloc
C
      REAL gmetr(4)
      SAVE gmetr
      INTEGER ifirst
      SAVE    ifirst
      DATA gmetr/-1.,-1.,-1.,1./
      DATA    ifirst/0/
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
      twopi=ACOS(0.)*4.
      itransl=0
C
      NEVSTAT(1)=NEVSTAT(1)+1
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
      amtarini=AMIN(2)
C
C---    The beam energy and the target 4-momenta
C
C---        No Fermi motion
C
      IF(NUCLTARG.EQ.0) THEN
         IF(IDBEAM.EQ.0.OR..NOT.HEXIST(IDBEAM+1)) GO TO 999
C
         IF(IBEATYP.EQ.0) THEN
            pbeam=EPH_LIM(1)    !   fixed energy
         ELSE
            ntry=0
 100        ntry=ntry+1
            pbeam=HRNDM1(IDBEAM+1)   ! beam spectrum * cross section
            IF(pbeam.LT.EPH_LIM(1).OR.pbeam.GT.EPH_LIM(2)) THEN ! check the boundaries
               IF(ntry.LT.1000) GO TO 100
               WRITE(lout,1000) pbeam,EPH_LIM
 1000          FORMAT(' *** ERROR: Beam energy out of range ',3F9.3)
               GO TO 999
            ENDIF
         ENDIF
         PIN(3,1)=pbeam
         NEVSTAT(6)=NEVSTAT(6)+1
      ELSE     !           Fermi motion (spectral functions)
         IF(IDBEAM.EQ.0.OR..NOT.HEXIST(IDBEAM+12)) GO TO 999
         IF(IDNUCLT.EQ.0.OR..NOT.HEXIST(IDNUCLT+11)
     +                  .OR..NOT.HEXIST(IDNUCLT+12)) GO TO 999
         IF(ifirst.EQ.0) THEN
            CALL HGIVE(IDNUCLT+12,tit,np,p1,p2,ne,e1,e2,nwt,lloc)
            ifirst=1
         ENDIF
C---       Simulate the momenta of the beam the target nucleon and randomly select the 
C---             combination duing the cross section energy dependence
         ntry0=0
 200     ntry0=ntry0+1
         DO i=1,2
            DO j=1,4
               pilab(j,i)=0.
            ENDDO
         ENDDO
         IF(IBEATYP.EQ.0) THEN
            pbeam=EPH_LIM(1)    !  beam  fixed energy
         ELSE
            pbeam=HRNDM1(IDBEAM+12)  ! beam random energy
         ENDIF
         pilab(3,1)=pbeam
         pilab(4,1)=SQRT(pbeam**2+AMIN(1)**2)
 210     CALL HRNDM2(IDNUCLT+11,psf,esf) ! coarse grid
         IF(psf.LT.p2.AND.esf.LT.e2) CALL HRNDM2(IDNUCLT+12,psf,esf) ! fine grid
         NEVSTAT(6)=NEVSTAT(6)+1
         psf=psf*0.001*SFSCALE   ! Mev --> GeV and scaling for tests
         esf=esf*0.001*SFSCALE   ! Mev --> GeV and scaling for tests
         eneff=AMNUCLS-SQRT((AMNUCLS-amtarini+esf)**2+psf**2)   ! energy of the struck nucleon
         ameff=eneff**2-psf**2
         if(ameff.LE.0.01) THEN    ! target mass >100 MeV - an arbitrary cut to avoid too large boosts to the target rest frame
            IF(ntry0.LT.2000) GO TO 210
            WRITE(lout,*) ' *** ERROR bg_eve low xsec, tries= ',ntry0
            GO TO 999
         ENDIF
         ameff=SQRT(ameff)
C
C---      PYTHIA only uses on-shell initial particles. We have to assume the struck proton is on shell
C
         IF(ISIMUL.EQ.0) THEN 
            ameff=AMIN(2)
            eneff=SQRT(ameff**2+psf**2)
         ENDIF
C
         pilab(4,2)=eneff
         ct=-1.+2.*RNDM(ct)
         st=SQRT(1.-ct**2)
         phi=twopi*RNDM(st)
         pilab(1,2)=psf*st*COS(phi)
         pilab(2,2)=psf*st*SIN(phi)
         pilab(3,2)=psf*ct
         ww=0.
         DO i=1,4
            ww=ww+gmetr(i)*(pilab(i,1)+pilab(i,2))**2
         ENDDO
         IF(ww.GT.0.) ww=SQRT(ww)
         xsec=HX(IDBEAM+13,ww)  ! cross section for this W (max element=1)
         xsrnd=RNDM(ww)
C         write(6,*) ' ===== ',pbeam,psf,esf
C         write(6,FMT='(4F12.4)') pilab
C         write(6,*) ww,xsec,ntry0
         IF(xsrnd.GT.xsec) THEN 
            IF(ntry0.LT.2000) GO TO 200
            WRITE(lout,*) ' *** ERROR bg_eve low xsec, tries= ',ntry0
            GO TO 999
         ENDIF
         AMIN(2)=ameff
C
C---       Go to the frame with the struck nucleon at rest 
C
         DO i=1,3
            bettar(i)=pilab(i,2)/pilab(4,2)
         ENDDO
         bettar(4)=pilab(4,2)/ameff
         CALL GLOREN(bettar,pilab(1,1),pitar(1,1))
C         CALL GLOREN(bettar,pilab(1,2),pitar(1,2))
         DO i=1,3
            pitar(i,2)=0.
            ptmp1(i,1)=pitar(i,1)   ! auxil - direction of beam in XY plane
         ENDDO
         pitar(4,2)=ameff
         ptmp1(3,1)=0.  ! auxil - direction of beam in XY plane
         
         CALL OMROTS1(pitar(1,1),ptmp1(1,1),rotz)  ! calculate rotation matrix, to beam alongz
         CALL OMROTV1(1,pitar(1,1),rotz,pitarr(1,1))
         pitarr(4,1)=pitar(4,1)
         DO i=1,4
            pitarr(i,2)=pitar(i,2)
         ENDDO
         itransl=1
C
C         write(6,FMT='('' pilab  '',4F12.4)') pilab
C         write(6,FMT='('' bettar '',4F12.4)') bettar
C         write(6,FMT='('' pitar  '',4F12.4)') pitar
C         write(6,FMT='('' pitarr '',4F12.4)') pitarr
C
C         DO i=1,3
C            bettar(i)=-bettar(i)
C         ENDDO
C         CALL OMROTV1(0,pitarr(1,1),rotz,pitar(1,1))
C         pitar(4,1)=pitarr(4,1)
C         DO i=1,4
C            pitar(i,2)=pitarr(i,2)
C         ENDDO
C         DO i=1,2
C            CALL GLOREN(bettar,pitar(1,i),pilab(1,i))
C         ENDDO
C         write(6,FMT='('' pitarr '',4F12.4)') pitarr
C         write(6,FMT='('' pitar  '',4F12.4)') pitar
C         write(6,FMT='('' pilab  '',4F12.4)') pilab
C
         DO i=1,2
            DO j=1,3
               PIN(j,i)=pitarr(j,i)
            ENDDO
         ENDDO
      ENDIF
C
      NTRA=0
C
C---   Choose the package
C
      ierr=0
      IF(ISIMUL.EQ.0) THEN
         IF(pbeam.LT.EPYMIN) THEN
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
C---     Rotate and boost into the Lab frame for the "Fermi" case
C
      IF(itransl.NE.0) THEN
         DO i=1,3
            bettar(i)=-bettar(i)
         ENDDO
         DO i=1,2
            DO j=1,3
               PIN(j,i)=pilab(j,i)
            ENDDO
         ENDDO
         DO ip=1,NTRA
            ww=AMTRA(ip)**2
            DO i=1,3
               potarr(i)=PTRA(i,ip)
               ww=ww+potarr(i)**2
            ENDDO
            potarr(4)=SQRT(ww)
            CALL OMROTV1(0,potarr(1),rotz,potar(1))
            potar(4)=potarr(4)
            CALL GLOREN(bettar,potar(1),polab(1))
            DO i=1,3
               PTRA(i,ip)=polab(i)
            ENDDO
         ENDDO
      ENDIF
C
      NEVSTAT(2)=NEVSTAT(2)+1
C
C---    Event selection (if needed)
C
      ierr=0
      IF(NFEFINAL.GT.0) THEN
         npfinal=0
         DO ip=1,NTRA
            IF(ITPTRA(2,ip).EQ.1) npfinal=npfinal+1 ! number of undecayed (final) particles
         ENDDO
         IF(npfinal.NE.NFEFINAL) THEN
            ierr=1
            GO TO 999
         ENDIF
      ENDIF
C
      DO i=1,MXEVSEL
         nkfsel=0
         IF(KFEVSEL(i).NE.0) THEN
            DO ip=1,NTRA
              IF(ITPTRA(3,ip).EQ.KFEVSEL(i)) THEN
                pt2=PTRA(1,ip)**2+PTRA(2,ip)**2
                ptot=SQRT(pt2+PTRA(3,ip)**2)
                IF(ptot.GT.0.01) THEN
                  thdg=ASIN(SQRT(pt2)/ptot)*180./3.1415
                  IF(ptot.GE.PAREVSEL(1).AND.ptot.LE.PAREVSEL(2).AND.
     +               thdg.GE.PAREVSEL(3).AND.thdg.LE.PAREVSEL(4)) THEN
                     nkfsel=nkfsel+1
                  ENDIF
                ENDIF
              ENDIF 
            ENDDO
            IF(nkfsel.EQ.0.OR.
     +        (NFEVSEL(i).GT.0.AND.nkfsel.NE.NFEVSEL(i))) THEN
               ierr=1           ! required particle not found in the record
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
C      write(6,*) pbeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C

