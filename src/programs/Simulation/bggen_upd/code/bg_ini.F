      SUBROUTINE BG_INI(IERR)
C
C---      Initialize the FFREAD and the relevant variables
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_reac.inc'
      INCLUDE 'bg_statc.inc'
      INCLUDE 'bg_tmp.inc'
      INCLUDE 'bg_beamext.inc'
C
      INTEGER mxffr,jffr
      PARAMETER (mxffr=10000)
      COMMON/CFREAD/ jffr(mxffr)
C
      INTEGER ier
C
      INTEGER i,j,lun,lout,iost,ip,kd(4),kf,lenc,idgea,ig,idec
     +       ,nw,iw(12),lw(12),kw(12)
      REAL am,wg,de
      CHARACTER cline*132,cnam*16
      INTEGER k_par(16)
      REAL    r_par(16)
      EQUIVALENCE (k_par(1),r_par(1))
      INTEGER nucleon     ! =0 - proton target, =1 - neutron target (for nuclei Spectral Function)
C
      lout=6
      IERR=1
C
      CALL HBOOK_INI
C
      NEVENT=0
      RUNNO=2
      IWROUT(1)=1
      IWROUT(2)=0
      IWROUT(3)=0
      IRND_SEQ=0
      NPRIEV=0

C
C---     Beam
C
      EPH_LIM(1)=0.15
      EPH_LIM(2)=12.
      EELEC=12.  ! electron energy
      EPEAK=9.   ! peak right edge
      EEMIT=1e-8 ! electron beam emittance (m.rad)
      RADT=50.E-6 ! radiator thickness (m)
      ZCOLL=7600.    ! Z distance to the collimator (cm)
      DCOLL=0.0050   ! collimator diameter (m)
      EPYMIN=3.      ! min energy for PYTHIA
      ELOWMIN=0.15   ! min energy for the low energy generator (GeV)
C
      IDBEAM=0
      NHBEA=0
      IBEATYP=-1
      DO i=1,2
         PHLIM(i)=0.
      ENDDO
C
      NBMEXT=0       ! if>0 points/bins in an externaly set beam energy distribution
      IBMEXT=0
      DO i=1,2
         EXT_LIM(i)=0.
      ENDDO
      DO i=1,MXBMPN
         RBEAEXT(i)=0.
         PBEAEXT(i)=0.
      ENDDO
C
C---               Standalone exclusive reactions
C
      ISIMUL=0       ! =0 - regular BG, >0 - exclusive reactions (=1 - J/psi)
      IPINIT(1)=1    ! beam particle type (GEANT) - photon 
      IPINIT(2)=14   ! target particle tyoe (GEANT) - proton
C
      LISTPYTH=0 ! .NE.0 print PYTHIA particle list at initialization
      AMLIST=6.  ! for LISTPYTH>0 print only particles with mass < AMLIST
C
      DO i=1,2
         IPREAC(i)=-1
      ENDDO
      ELREAC(1)=-1.
      ELREAC(2)=-1.
      NPXREAC=0
      DO i=1,MXPNTR
         XSREAC(i)=0.
         ESREAC(i)=0.
      ENDDO
      TLREAC(1)=0.
      TLREAC(2)=0.
      DO i=1,MXTRNG
         DO j=1,11
            TPAREAC(j,i)=0.
         ENDDO
         ITDRNG(i)=0
      ENDDO
C
C---      Particle selection
C
      DO i=1,MXEVSEL
         KFEVSEL(i)=0      ! types
         NFEVSEL(i)=0      ! types
      ENDDO
      NFEFINAL=0           ! # of particles in the final state KS=1
      PAREVSEL(1)=0.       ! Pmin
      PAREVSEL(2)=99999.   ! Pmax
      PAREVSEL(3)=0.       ! theta_min
      PAREVSEL(4)=180.     ! Theta_max
C
      DO i=1,MXSTAT
         NEVSTAT(i)=0
      ENDDO
C
      DO i=1,MXRADP
         IRADCORR(i)=0             ! Rad corrections for decays (=0 - off)
      ENDDO
C !!!      IRADINI=0
C
      LUNWR(1)=0   ! HDDS file - LUN not used
      LUNWR(2)=2   ! sequential file
      LUNWR(3)=3   ! NTUPLE file
C
C
C---   Read the particle masses (GEANT numbering)
C---   Read the GEANT<->PYTHIA particle table
C
      NPADEF=0
      DO ig=1,MXPART
         KF_PART(ig)=0
         KC_PART(ig)=0
         LD_PART(ig)=0
         AM_PART(ig)=0.
         WG_PART(ig)=0.
         ND_PART(ig)=0
         DO i=1,4
            KD_PART(i,ig)=0
         ENDDO
         CH_PART(ip)='#'
      ENDDO
C
      DO i=1,20          !  reset the temp. arrays
         K_PTMP(i)=0
      ENDDO
      DO i=1,4
         CALL UCTOH('    ',K_PTMP(9+i),4,4)
      ENDDO
C
      DO i=-MXPKC,MXPKC
         IDG_KC(i)=0
         KF_KC(i)=0
         MULT_KC(i)=0
         NODC_KC(i)=0
C         CH_KC(i)=' '
      ENDDO
C
      NPYMOD=0
      DO i=1,8
         CALL UCTOH('    ',KPYTMP(i),4,4)
      ENDDO
C
      NUCLTARG=0     ! nuclear target A for quasielastic processes
      IGNUCLT=0
      SFSCALE=1.     ! scale factor for the spectral distribution
C
C---    Redefine FFREAD settings
C
      CALL FFINIT(mxffr)
      CALL FFSET('LINP',15)
      CALL FFSET('LOUT',6)
      CALL FFSET('SIZE',16)
      CALL FFSET('LENG',120)
C
C---    BGGEN keys
C
      CALL FFKEY('TRIG'      , NEVENT     , 1,'INTEGER')
      CALL FFKEY('RUNNO'     , RUNNO      , 1,'INTEGER')
      CALL FFKEY('WROUT'     , IWROUT(1)  , 3,'INTEGER')
      CALL FFKEY('RNDMSEQ'   , IRND_SEQ   , 1,'INTEGER')
      CALL FFKEY('NPRIEV'    , NPRIEV     , 1,'INTEGER')
      CALL FFKEY('LISTPYTH'  , LISTPYTH   , 1,'INTEGER')
      CALL FFKEY('LISTMLIM'  , AMLIST     , 1,'REAL')
      CALL FFKEY('EPHLIM'    , EPH_LIM(1) , 2,'REAL')
      CALL FFKEY('EELEC'     , EELEC      , 1,'REAL')
      CALL FFKEY('EPEAK'     , EPEAK      , 1,'REAL')
      CALL FFKEY('EEMIT'     , EEMIT      , 1,'REAL')
      CALL FFKEY('RADTHICK'  , RADT       , 1,'REAL')
      CALL FFKEY('ZCOLLIM'   , ZCOLL      , 1,'REAL')
      CALL FFKEY('DCOLLIM'   , DCOLL      , 1,'REAL')      
      CALL FFKEY('EPYTHMIN'  , EPYMIN     , 1,'REAL')
      CALL FFKEY('ELOWMIN'   , ELOWMIN    , 1,'REAL')
      CALL FFKEY('VERTEX'    , VERTEX(1)  , 3,'REAL')
      CALL FFKEY('SIMUL'     , ISIMUL     , 1,'INTEGER')
      CALL FFKEY('PARTINIT'  , IPINIT(1)  , 2,'INTEGER')
      CALL FFKEY('RADCORR'   , IRADCORR(1), MXRADP,'INTEGER')
      CALL FFKEY('REACPART'  , IPREAC(1)  , 2,'INTEGER')
      CALL FFKEY('REACTLIM'  , TLREAC(1)  , 2,'REAL')
      DO i=1,MXTRNG
         WRITE(cnam,FMT='(''REACTDIS'',I1.1)') i
         CALL FFKEY(cnam(1:9), TPAREAC(1,i),11,'REAL')
      ENDDO
      CALL FFKEY('REACENER'  , ELREAC(1)  , MXPNTR+3,'MIXED')
      CALL FFKEY('SELECTEV'  , KFEVSEL(1) , MXEVSEL,'INTEGER')
      CALL FFKEY('SELECTEVN' , NFEVSEL(1) , MXEVSEL,'INTEGER')
      CALL FFKEY('SELECTKIN' , PAREVSEL(1),  4,'REAL')
      CALL FFKEY('NPARFINAL' , NFEFINAL   ,  1,'INTEGER')
      CALL FFKEY('PARTICLE'  , K_PTMP(1)  , 16,'MIXED')
      CALL FFKEY('PYTMOD'    , KPYTMP(1)  ,  8,'MIXED')
      CALL FFKEY('NBEAMEXT'  , NBMEXT     ,  4,'MIXED')
      CALL FFKEY('RBEAMEXT'  , RBEAEXT(1) , MXBMPN,'REAL')
      CALL FFKEY('PBEAMEXT'  , PBEAEXT(1) , MXBMPN,'REAL')
      CALL FFKEY('TNUCLEUS'  , NUCLTARG   ,  1,'INTEGER')
      CALL FFKEY('TNUSCALE'  , SFSCALE    ,  1,'REAL')
C
      CALL FFGO
C
      ier=0
      DO ig=1,MXPART
         IF(KF_PART(ig).NE.0) THEN
            DO i=1,3
               IF(KD_PART(i,ig).GT.0) THEN
                  IF(KD_PART(i,ig).LE.MXPART) THEN
                     IF(KF_PART(KD_PART(i,ig)).NE.0) THEN
                        ND_PART(ig)=ND_PART(ig)+1
                     ELSE
                        ier=ig
                     ENDIF   
                  ELSE
                     ier=ig
                  ENDIF
               ELSE
                  GO TO 14
               ENDIF
            ENDDO
 14         IF(ND_PART(ig).GT.0.AND.ND_PART(ig).LT.2) ier=ig
C            WRITE(6,*) ig,AM_PART(ig),ND_PART(ig),KD_PART(1,ig)            
         ENDIF
      ENDDO
      IF(ier.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: in particle data, ID=',ier
         GO TO 999
      ENDIF
C
      DO i=1,2
         ier=1
         IF(IPINIT(i).GT.0.AND.IPINIT(i).LE.MXPART) THEN
            IF(KF_PART(IPINIT(i)).NE.0) THEN
               ier=0
            ENDIF
         ENDIF
         IF(ier.NE.0) THEN
            WRITE(6,1001) i,IPINIT(i)
 1001              FORMAT('  *** Init. error: PARTINIT:',I1,I5,
     +                                    ' not defined')
            GO TO 999
         ENDIF
      ENDDO
C
      IF(NUCLTARG.GT.0) THEN  ! find the target nuclei (optional)
         DO ig=1,MXPART
            IF(KF_PART(ig).GT.1000000000) THEN   ! nuclei definitions
               i=MOD(KF_PART(ig),10000)/10  ! KF code is 100ZZZAAAI
               IF(i.EQ.NUCLTARG) IGNUCLT=ig
            ENDIF
         ENDDO
         IF(IGNUCLT.GT.0) THEN
            AMNUCLS=AM_PART(IGNUCLT)
            AMNUCLN=AMNUCLS/NUCLTARG
         ELSE
            AMNUCLN=0.931  ! average nucleon mass inside a nucleus
            AMNUCLS=AMNUCLN*NUCLTARG
         ENDIF
      ENDIF
C
C---      Print particle table
C
      WRITE(lout,1000) 
 1000 FORMAT(1X,80('='),/,' Particle list locally defined '
     +      ,'  ID_GEant<33 are also defined in the GEANT3 code'
     +      ,/,' ID_Geant    KF PY_dec N_dec  ID_decays'
     +      ,' dec_flg     mass     width     name')
      DO ig=1,MXPART
         IF(KF_PART(ig).NE.0) THEN
            WRITE(lout,1002) ig,KF_PART(ig),LD_PART(ig),ND_PART(ig)
     +                ,(KD_PART(i,ig),i=1,4),AM_PART(ig),WG_PART(ig)
     +                ,CH_PART(ig) 
 1002       FORMAT(1X,I3,1X,I10,3X,I3,3X,I2,2X,3I4,3X,I2,3X,F10.6,2X
     +            ,F7.4,5X,A16)
         ENDIF
      ENDDO
      WRITE(lout,FMT='(1X,80(''-''))')
C
      CALL RND_INI(IRND_SEQ)  ! random number initialization
C
      IDBEAM=9000
      CALL BEAM_INI(ier)      ! initialize the beam flux (to be multiplied by the cross section later)
      IF(ier.NE.0) GO TO 999
C
      IFPYTH=0   ! PYTHIA flag
      IDLOWEN=0  ! low energy BG reactions flag
      IDTDIST=0  ! Excl. reaction ID for hist
C
C---     A nuclear target: read the spectral function
C
      IF(NUCLTARG.GT.0) THEN
         IDNUCLT=11000
         kf=KF_PART(IPINIT(2))  
         IF(kf.EQ.2212) THEN
            nucleon=0   ! proton
         ELSE IF(kf.EQ.2112) THEN
            nucleon=1   ! neutron
         ELSE
            ier=1
            WRITE(lout,1010) NUCLTARG,IPINIT(2)
 1010       FORMAT(' *** BG_INIT error: nucleus target A=',I4
     +          ,' interation on ig=',I4,' not proton or neutron')
            GO TO 999
         ENDIF
         CALL SF_INIT(NUCLTARG,IDNUCLT,nucleon,ier)
         IF(ier.NE.0) GO TO 999
      ENDIF
      IF(NUCLTARG.GT.0) THEN
         WRITE(lout,*) ' --- Target nucleus A=',NUCLTARG
     +                ,' GEANT type=',IGNUCLT,' ID=',IDNUCLT
      ENDIF
C
C---    Check the type of simulation - PYTHIA versus exclusive reactions
C
      IF(ISIMUL.EQ.0) THEN  ! BG simulation
C
C---   Pythia
C
         IF(EPH_LIM(2).GT.EPYMIN) THEN
            CALL PYTH_INI(ier)
            IF(ier.NE.0) GO TO 999
            IFPYTH=1
         ENDIF
C
C---   Low energy processes
C
         IF(EPH_LIM(1).LT.EPYMIN) THEN
            IDLOWEN=10000
            CALL LOWEN_INI(ier)
            IF(ier.NE.0) GO TO 999
         ENDIF
C
      ELSE
C
C---      Exclusive reactions
C
         IDTDIST=8000
         CALL REAC_INI(ier)
         IF(ier.NE.0) GO TO 999
      ENDIF
C
C---    The "effective" spectrum - the beam spectrum times the cross section
C
      CALL BEAM_EFF(ier)
      IF(ier.NE.0) GO TO 999
C
C---      Rad corrections
C
      NRADCOR=0
      DO i=1,MXRADP
         idgea=IRADCORR(i)
         IF(idgea.LE.0.OR.idgea.GT.MXPART) GO TO 100
         NRADCOR=NRADCOR+1
      ENDDO
 100  CONTINUE
      IF(NRADCOR.GT.0) THEN
         CALL PHOINI            ! Initialize PHOTOS for radiative corrections
      ENDIF
C
C---   Output file for HDDM
C
      IF(IWROUT(1).NE.0) THEN
         CALL OPEN_HDDM_OUTPUT('bggen.hddm')
      ENDIF
C
C---   Sequential output file
C
      IF(IWROUT(2).NE.0) THEN
         OPEN(LUNWR(2),FILE='bggen.dat',STATUS='UNKNOWN'
     +       ,FORM='UNFORMATTED')
      ENDIF
C
C---   NTUPLE
C
      IF(IWROUT(3).NE.0) THEN
         CALL BG_NTUP_INI(ier)
         IF(ier.NE.0) GO TO 999
      ENDIF
C
      flush(6)
      IERR=0
 999  RETURN
      END

      
