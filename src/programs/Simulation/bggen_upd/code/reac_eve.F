      SUBROUTINE REAC_EVE(IERR)
C
C---  Simulates 1 event - a single reaction
C     Reaction: gamma+p --> A+B , A - recoil, B - meson 
C               ISIMUL=1 A=p (14), B=J/psi (83)
C
      IMPLICIT NONE
      INTEGER IERR
C
      INCLUDE 'bg_ctrl.inc'
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_evec.inc'
      INCLUDE 'bg_reac.inc'
      INCLUDE 'hepevt.inc'
C
      REAL RNDM,GBRWIGN,HRNDM1
C
      INTEGER i,j,ip,np,idgea,ntry
     +       ,idgd,ntry1,ihel,ir
      REAL ebeam,ecm,ecm2,bet(4),qq,ct,st,phi
     +    ,twopi
     +    ,amtot    ! sum of the masses 
     +    ,ppf,epf1,ppi,tt,tmn,tmx,amdec
     +    ,wdm
C
      CHARACTER cent(2)*6,cproc*16
C
      INTEGER mxoutl
      PARAMETER (mxoutl=50)
      REAL ami(2),pcmi(4,2),plabi(4,2)
     +     ,am(mxoutl),pcm(4,mxoutl),plab(4,mxoutl),amd(6,mxoutl)
      INTEGER idg(mxoutl)        ! Geant ID
     +        ,iamf(mxoutl)      ! >0 if the mass has been defined (to address resonances)
     +        ,itorig(mxoutl)
     +        ,it1dec(mxoutl)
C
C     ------------------------------------------------------------------
C
C      tslreac=1.2  ! test, commented out
      IERR=1
      NTRA=0
      IF(ISIMUL.LT.1.OR.ISIMUL.GT.1) GO TO 999
C
C---    Beam energy
C
      ebeam=PIN(3,1)
C
C---      Initial state (beam goes along Z - no rotation applied)
C
      DO i=1,2
         ami(i)=AMIN(i)
         qq=0.
         DO j=1,3
            plabi(j,i)=PIN(j,i)
            qq=qq+plabi(j,i)**2
         ENDDO
         plabi(4,i)=SQRT(qq+ami(i)**2)
      ENDDO
C
C      write(6,*) 'ami', ami,plab(4,1),plab(4,2)
      qq=plabi(4,1)*plabi(4,2)
      DO i=1,3
         qq=qq-plabi(i,1)*plabi(i,2)
      ENDDO
      ecm2=ami(1)**2+ami(2)**2+2.*qq
      ecm=SQRT(ecm2)
C
      DO i=1,3
         bet(i)=(plabi(i,1)+plabi(i,2))/(plabi(4,1)+plabi(4,2))
      ENDDO
      bet(4)=(plabi(4,1)+plabi(4,2))/ecm
C
      DO i=1,2
         CALL GLOREN(bet,plabi(1,i),pcmi(1,i))
      ENDDO
      DO i=1,3
         bet(i)=-bet(i)
      ENDDO
C
      DO i=1,mxoutl
         idg(i)=0
         iamf(i)=0
      ENDDO
C
      ntry=0
 30   ntry=ntry+1
      amtot=0.
      DO ip=1,2
         idgea=IPREAC(ip)  ! defined, already checked 
         idg(ip)=idgea
         itorig(ip)=0
         it1dec(ip)=0
         iamf(ip)=0
         ntry1=0
 35      ntry1=ntry1+1       ! find the combination of masses (mother and daughters) - energy conserv.
         am(ip)=AM_PART(idgea)
         IF(WG_PART(idgea).GT.0.) 
     +                   am(ip)=am(ip)+WG_PART(idgea)*GBRWIGN(am)         
         amdec=0.
         IF(ND_PART(idgea).GT.1) THEN 
            DO i=1,ND_PART(idgea)
               idgd=KD_PART(i,idgea) ! defined, already checked
               amd(i,ip)=AM_PART(idgd)
               IF(WG_PART(idgd).GT.0.) 
     +                    amd(i,ip)=amd(i,ip)+WG_PART(idgd)*GBRWIGN(amd)         
               amdec=amdec+amd(i,ip)
            ENDDO
         ENDIF
         IF(am(ip).LT.amdec+0.01) THEN
            IF(ntry1.LT.1000) GO TO 35
            WRITE(6,*) ' *** REAC_EVE unsuff mass for decay '
     +           ,idgea,am(ip),amdec 
            GO TO 999
         ENDIF
         iamf(ip)=1
         amtot=amtot+am(ip)
      ENDDO
C
      IF(amtot.GE.ecm-0.001) THEN
         IF(ntry.LT.1000) GO TO 30
         WRITE(6,*) ' *** REAC_EVE unsuff energy for reaction '
     +           ,amtot,ecm 
         GO TO 999
      ENDIF
      np=2
C
C---  Treat the kinematics as 2-body one, in CM
C
      twopi=ACOS(0.)*4.
C
C---      In CM: momentum and energies of the particles
C
      epf1=(ecm2+am(1)**2-am(2)**2)/2./ecm
      qq=epf1**2-am(1)**2
      IF(qq.LE.0.) GO TO 999
      ppf =SQRT(qq)             ! final momentum in CM
      ppi =SQRT(pcmi(4,2)**2-ami(2)**2) ! initial momentum in CM
C
      qq=ami(1)**2+am(1)**2-2.*epf1*pcmi(4,1)
      tmx=qq+2.*ppf*ppi
      tmn=qq-2.*ppf*ppi
      IF(ITDSING.NE.0) THEN
         IF(TSLREAC.LT.0.001) THEN
            tt=tmn+(tmx-tmn)*RNDM(qq)
         ELSE
            tt=1./TSLREAC*ALOG(EXP(TSLREAC*tmn)
     +        +RNDM(qq)*(EXP(TSLREAC*tmx)-EXP(TSLREAC*tmn)))  ! random tt - a single exponential
         ENDIF
      ELSE      ! Find the proper energy range for the t-dep functions
         DO ir=1,MXTRNG
            IF(ITDRNG(ir).GT.0) THEN
               IF(ebeam.LE.TPAREAC(1,ir)) THEN
                  ntry=0
 40               ntry=ntry+1
                  tt=HRNDM1(IDTDIST+10*ir) ! random tt - a stored function
                  IF(tt.GE.tmn.AND.tt.LE.tmx) THEN
                     GO TO 50
                  ELSE
                     IF(ntry.LT.1000) GO TO 40
                     WRITE(6,*) ' *** REAC_EVE t out of range',tt
                     GO TO 999
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
 50      CONTINUE
      ENDIF
      ct=(tt-qq)/2./ppf/ppi
C
      st=SQRT(1.-ct**2)
      phi=twopi*RNDM(st)
C
C---    2-body
C
      pcm(4,1)=epf1
C     
      pcm(1,1)=ppf*st*COS(phi)
      pcm(2,1)=ppf*st*SIN(phi)
      pcm(3,1)=ppf*ct
C
      DO i=1,3
         pcm(i,2)=-pcm(i,1)
      ENDDO
      pcm(4,2)=ecm-pcm(4,1)
C
C---       Boost to Lab
C
      DO i=1,2
         CALL GLOREN(bet,pcm(1,i),plab(1,i))
      ENDDO
C
C---       Decays?
C
      ip=0
 100  ip=ip+1     ! check the whole list of generated particles: a loop
      idgea=idg(ip)
      IF(ND_PART(idgea).GT.0) THEN
C         write(6,*) 'ip,idgea,ndec ',ip,idgea,ND_PART(idgea)
         it1dec(ip)=np+1
         ntry=0
 120     ntry=ntry+1
         amdec=0.
         DO i=1,ND_PART(idgea)
            idgd=KD_PART(i,idgea)
            idg(np+i)=idgd
            IF(iamf(ip).EQ.0) THEN
               amd(i,ip)=AM_PART(idgd)  
               IF(WG_PART(idgd).GT.0.) 
     +                    amd(i,ip)=amd(i,ip)+WG_PART(idgd)*GBRWIGN(amd)         
            ENDIF
            amdec=amdec+amd(i,ip)            
            am (np+i)=amd(i,ip) 
            itorig(np+i)=ip
            it1dec(np+i)=0
         ENDDO
         IF(am(ip).LT.amdec+0.01) THEN
            IF(ntry.LT.1000) GO TO 120
            WRITE(6,*) ' *** REAC_EVE unsuff mass for decay '
     +           ,idgea,am(ip),amdec 
            GO TO 999
         ENDIF
         iamf(ip)=1
         IF(ND_PART(idgea).EQ.2) THEN       ! 2-body decay
            ihel=KD_PART(4,idgea)     ! decay angle flag =0 - unoform, =1 - rho-like, =2 - j/psi-like
            CALL OMDECA2(plab(1,ip),amd(1,ip),ihel,plab(1,np+1))
         ELSE IF(ND_PART(idgea).EQ.3) THEN
            CALL OMDECA3(plab(1,ip),amd(1,ip),0.,plab(1,np+1))
         ENDIF
         np=np+ND_PART(idgea)
      ENDIF
      IF(ip.LT.np) GO TO 100   ! Loop for all entries
C
      DO ip=1,np
         idgea=idg(ip)
         DO j=1,3
            PTRA(j,ip)=plab(j,ip)
         ENDDO
         AMTRA(ip)=am(ip)
         ITPTRA(1,ip)=idgea
         DO j=2,6
            ITPTRA(j,ip)=0
         ENDDO
C         write(6,*) i,idg(ip),MXPART,KF_PART(idg(ip)),itorig(ip),it1dec(ip)
         IF(idgea.GT.0.AND.idgea.LE.MXPART) THEN
            ITPTRA(3,ip)=KF_PART(idgea)
         ENDIF
         ITPTRA(4,ip)=itorig(ip)
         ITPTRA(5,ip)=it1dec(ip)
         IF(it1dec(ip).GT.0) ITPTRA(6,ip)=it1dec(ip)+ND_PART(idgea)-1
         ITPTRA(2,ip)=1
         IF(it1dec(ip).NE.0) ITPTRA(2,ip)=10   ! indicates that this particle should not be used in GEANT 
      ENDDO
      NTRA=np
C
      IERR=0
 999  CONTINUE
C      write(6,*) ebeam,IEVPROC,ibin,xstot,xssum,NTRA
C
      END
C

