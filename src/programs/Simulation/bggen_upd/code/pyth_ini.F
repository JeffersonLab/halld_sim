      SUBROUTINE PYTH_INI(IERR)
C
C---  Initilize PYTHIA
C     Reaction: gamma+p
C---  Input: /phctrl/ beam energy 
C            file "pythia-geant.dat" contains a table for PYTHIA<->GEANT particle code conversion
C            file "pythia.dat" - redefinition of PYTHIA parameters (from HERMES, adapted to GLUEX)           
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_tmp.inc'
C
      INTEGER PYCOMP
      EXTERNAL PYCOMP
C
      DOUBLE PRECISION dbeam
      CHARACTER cpar*100,cpini(2)*16
      REAL am(2),ss
C
      INTEGER lun,lout,i,j,lenc,kf,kc,ks,iost,ig
C
C     ------------------------------------------------------------------
C
      IERR=1
C
      lun=9
      lout=6
C
C---    KF/KC/GEANT mapping
C
      DO ig=1,MXPART
         kf=KF_PART(ig)
         IF(kf.NE.0) THEN
            kc=PYCOMP(kf)
            IF(kc.GT.0.AND.kc.LE.MXPKC) THEN
               IF(kf.LT.0) kc=-kc
               KC_PART(ig)=kc
               IDG_KC(kc)=ig
               kc=ABS(kc)
C                  write(lout,FMT='(10I8)') ig,kf,kc
C
C---             Forbid the decays for particles if the proper flag set
C
               IF(LD_PART(ig).EQ.0) THEN
                  WRITE(cpar,1000) kc,0
 1000             FORMAT('MDCY(',I4,',1)=',I2)
                  CALL PYGIVE(cpar)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
C 
C---   Read the pythia settings for JLab energies
C      
      OPEN(lun,FILE='pythia.dat',STATUS='OLD',IOSTAT=iost
     +       ,FORM='FORMATTED')
      IF(iost.NE.0) THEN
         WRITE(lout,*) ' *** ERROR: Missing file pythia.dat'
         GO TO 999
      ENDIF
 30   READ(lun,'(A)',IOSTAT=iost) cpar
      IF(iost.EQ.0) THEN
         CALL PYGIVE(cpar)
         GO TO 30
      ELSE IF(iost.GT.0) THEN
         WRITE(lout,*) ' *** ERROR: Reading file pythia.dat'
         GO TO 999
      ENDIF
      CLOSE(lun)
C
      DO i=1,NPYMOD
C         write(6,FMT='(I4,2X,A)') i,CPYMOD(i)
         CALL PYGIVE(CPYMOD(i)) ! PYTHIA modifications from FFREAD
      ENDDO
C
C---      Variable energy?
C
      IF(NHBEA.GT.1) THEN
         cpar='MSTP(171)=1'
         CALL PYGIVE(cpar)
         cpar='MSTP(172)=1'
         CALL PYGIVE(cpar)
      ENDIF
C
C---   Beam and target
C
      DO i=1,2
         ig=IPINIT(i)
         am(i)=AM_PART(ig)
         kf=KF_PART(ig)     ! IPINIT validity has been checked in bg_ini.F
         CALL PYNAME(kf,cpini(i))
      ENDDO
      ss=SQRT(am(1)**2+am(2)**2+2.*EPYMIN*am(2))  ! min sqrt(s); target at rest
      WRITE(cpar,FMT='(''PARP(2)='',F4.2)') ss*0.95
      CALL PYGIVE(cpar)       ! sets the minimal PYTHIA CM energy, with a 5% margin (without a margin the job may stop)  
C
C---   Initialize PYTHIA
C
      dbeam=DBLE(EPH_LIM(2))
C
      CALL PYINIT('FIXT',cpini(1),cpini(2),dbeam)
C
      IF(LISTPYTH.EQ.-11.OR.LISTPYTH.EQ.-12) CALL PYLIST(-LISTPYTH)  ! print a list of particles
      IF(LISTPYTH.GT.0) CALL PYTH_PRI(LISTPYTH,AMLIST)               ! print a list of particles
C
      IERR=0
 999  CONTINUE
C
      END

