C=======================================================================
      SUBROUTINE SF_INIT(NUCLA,IDSF,NUCLEON,IERR)
C      SUBROUTINE SF_INIT(NUCLEON,ID1,ID2,ANUCL,ierr)
C
C---    Initialize the Spectral Function generator: read a file sf_nucleus.dat
C---    INPUT:
C---      NUCLA     - A (integer)
C---      NUCLEON=0 - protons, =1 - neutrons
C---      IDSF - ID for histograms: IDSF+1 - coarse, IDSF+2 - fine; +11, +12 - for HRNDM
C---    OUTPUT:
C---      IERR !=0 - error
C
      IMPLICIT NONE
      INTEGER NUCLA,NUCLEON,IDSF,IERR
C
      INCLUDE 'sf_data.inc'
C
      CHARACTER cfile*72,cline*120
C
      LOGICAL HEXIST
      EXTERNAL HEXIST
      REAL SF_LOOKUP,HSUM
      EXTERNAL SF_LOOKUP
C
      INTEGER lun
      INTEGER iost
      INTEGER nlines        ! # of data lines read (without the header) 
     +       ,nltot         ! # total number of lines read
      INTEGER nuclaf        ! atomic mass A read from the sf file
     +       ,ie,ip,ne,np,npi,nei,i,i1
C
      REAL pp(6)
      DOUBLE PRECISION dnorm
C
      REAL p1,p2,e1,e2,p2i,e2i,sf,ppc,eec,pb,eb,pbi,ebi,sfsum
      CHARACTER tit*80,cnm*4,ca*3
      INTEGER nav,iave,iavp
      DOUBLE PRECISION dav
      DATA nav/10/
C
C     ------------------------------------------------------------------
C
      IERR=1
      lun=9
C
      IF(IDSF.LT.1) THEN
         WRITE(6,*) ' SF_INIT error IDSF= ',IDSF
         GO TO 999
      ENDIF
C
      IF(HEXIST(IDSF+1).OR.HEXIST(IDSF+2)) THEN
         WRITE(6,*) ' SF_INIT error IDSF+1 or IDSF+2 exist, IDSF= ',IDSF
         GO TO 999
      ENDIF         
C
      IF(NUCLA.LT.1.OR.NUCLA.GT.999) THEN
         WRITE(6,*) ' SF_READ: wrong A=',NUCLA
         GO TO 999
      ENDIF
C
      IF(NUCLA.LT.10) THEN
         i1=1
         WRITE(ca,FMT='(I1)') NUCLA
      ELSE IF(NUCLA.LT.100) THEN
         i1=2
         WRITE(ca,FMT='(I2)') NUCLA
      ELSE
         i1=3
         WRITE(ca,FMT='(I3)') NUCLA
      ENDIF
C
C      cfile='sf_nucleus.dat'
      cfile='spec_fun/sf_Benhar_'//ca(1:i1)//'.dat'
C      write(6,*) ' Open file ',cfile
      OPEN(UNIT=lun,FILE=cfile,STATUS='OLD',IOSTAT=iost)
      IF(iost.NE.0) THEN
         WRITE(6,*) ' SF_INIT: error opening ',cfile
         GO TO 999
      ENDIF
      REWIND lun
      WRITE(6,*) ' --- Nuclear target, spectral function file ',cfile
C
      NE_SF=0
      NP_SF=0
C
      nlines=0
      nltot=0
      nuclaf=0
      np=0
      ne=0
C
 100  READ(lun,FMT='(A)',ERR=997,END=998) cline          ! read the header
C        write(6,*) cline
        nltot=nltot+1
        IF(cline(1:1).NE.'*') GO TO 110
        WRITE(6,*) ' ---    Header ',cline(1:70)
        i=INDEX(cline,'A=') 
        IF(i.GT.1) READ(cline(i+2:i+6),*) nuclaf
      GO TO 100
C
 110  CONTINUE
      IF(nuclaf.NE.NUCLA) THEN
         WRITE(6,*) ' SF_INIT: error: NUCLA=',NUCLA
     +             ,' but the file belongs to A=',nuclaf
         GO TO 999
      ENDIF
C
      READ(cline,*) np,ne
      nlines=nlines+1
C      WRITE(6,*) ' --- Spectral function defined, ne,np = ',ne,np
C
      IF(ne.LT.1.OR.np.LT.1.OR.nuclaf.LT.1) THEN
         WRITE(6,*) ' SF_INIT: error reading header ',np,ne,nuclaf
         GO TO 999
      ENDIF
      IF(ne.GT.MXSFB.OR.np.GT.MXSFB) THEN
         WRITE(6,*) ' SF_INIT: error np,ne > MXSFB   ',np,ne,nuclaf
         GO TO 999
      ENDIF
C           
      dnorm=0.D0  ! normalization
C
      DO ip=1,np
        DO ie=1,ne
C          write(6,*) ip,ie,nlines
          READ(lun,*,ERR=997,END=998) pp
          nlines=nlines+1
          nltot=nltot+1
          IF(ip.EQ.1) THEN
             SF_ECEN(ie)=pp(2)
             SF_EWID(ie)=pp(6)
          ELSE
             IF(ABS(pp(2)-SF_ECEN(ie)).GT.0.0001) THEN
                WRITE(6,*) ' SF_INIT: error line,ip,ie,ec(ie)!=pp(2)'
     +                    ,nlines,ip,ie,SF_ECEN(ie),pp
                GO TO 999
             ENDIF
          ENDIF
          IF(ie.EQ.1) THEN
             SF_PCEN(ip)=pp(1)
             SF_PWID(ip)=pp(5)
          ELSE
             IF(ABS(pp(1)-SF_PCEN(ip)).GT.0.0001) THEN
                WRITE(6,*) ' SF_INIT: error line,ip,ie,pc(ip)!=pp(1)'
     +                    ,nlines,ip,ie,SF_PCEN(ip),pp
                GO TO 999
             ENDIF
          ENDIF
C
          IF(NUCLEON.EQ.0) THEN
             SF_PROB(ie,ip)=pp(3)   ! for protons
          ELSE
             SF_PROB(ie,ip)=pp(4)   ! for neutrons
          ENDIF
          IF(SF_PROB(ie,ip).LT.0.) THEN
             WRITE(6,*) ' SF_INIT: error line,ip,ie,sf<0 '
     +                ,nlines,ip,ie,SF_PROB(ie,ip)
             SF_PROB(ie,ip)=0.
          ENDIF
          dnorm=dnorm+DBLE(SF_PROB(ie,ip))
          SF_PROB(ie,ip)=SF_PROB(ie,ip)/pp(5)/pp(6)   ! diff. probab. dP(e,p)/dedp
C          IF(ip.EQ.10) WRITE(6,*) ip,ie,pp,SF_PROB(ie,ip)
C
        ENDDO
      ENDDO
      CLOSE(UNIT=lun)
C
C---   Normalization
C
      DO ip=1,np
        DO ie=1,ne
           SF_PROB(ie,ip)=SF_PROB(ie,ip)/REAL(dnorm)
C           WRITE(6,2010) SF_PCEN(ip),SF_ECEN(ie),SF_PROB(ie,ip)
C     +                  ,SF_PWID(ip),SF_EWID(ie)
C 2010      FORMAT(2F9.3,2X,E12.4,2X,2F9.3)
        ENDDO
      ENDDO
C
      NE_SF=ne
      NP_SF=np
      SF_A=nuclaf
      WRITE(6,*) ' --- Spectral function defined for A=',nuclaf
     +           ,'  ne,np= ',NE_SF,NP_SF
     +           ,' initial normalization=',REAL(dnorm)
C      DO ip=1,np
C        DO ie=1,ne
C           WRITE(6,2010) SF_PCEN(ip),SF_ECEN(ie),SF_PROB(ie,ip)
C     +                  ,SF_PWID(ip),SF_EWID(ie)
C 2010      FORMAT(2F9.3,2X,E12.4,2X,2F9.3)
C        ENDDO
C      ENDDO
C
C---   Book histograms: fill the probablility/bin in each bin: d^2(P(E,P))/de/dp * de*dp
C
      p1=0.
      pb=10.  ! 10 MeV bin   
      p2=SF_PCEN(NP_SF)+SF_PWID(NP_SF)/2.
      np=INT((p2-p1)/pb)+1   
      p2=p1+pb*np
C        
      e1=0.
      eb=10.  ! 10 MeV bin   
      e2=SF_ECEN(NE_SF)+SF_EWID(NE_SF)/2.
      ne=INT((e2-e1)/eb)+1   
      e2=e1+eb*ne
C
      pbi=2.  ! 2MeV - fine binning
      ebi=2.
      npi=200
      nei=100
      p2i=p1+npi*pbi
      e2i=e1+nei*ebi
C
      write(6,*) np,p2,pb,ne,e2,eb
      write(6,*) npi,p2i,pbi,nei,e2i,ebi
C
      WRITE(cnm,FMT='(I3)') nuclaf
      tit='Spectral function (coarse) for A='//cnm
      CALL HBOOK2(IDSF+1,tit,np,p1,p2,ne,e1,e2,0.)
      tit='Spectral function (fine) for A='//cnm
      CALL HBOOK2(IDSF+2,tit,npi,p1,p2i,nei,e1,e2i,0.)
C
      DO ip=1,np
         DO ie=1,ne
            dav=0.D0
            DO iavp=1,nav
               DO iave=1,nav
                  ppc=p1+(ip-1.+(iavp-0.5)/REAL(nav))*pb
                  eec=e1+(ie-1.+(iave-0.5)/REAL(nav))*eb
                  sf=SF_LOOKUP(ppc,eec)
                  dav=dav+DBLE(sf)
               ENDDO
            ENDDO
            sf=REAL(dav)/nav**2*pb*eb
            CALL HF2(IDSF+1,ppc,eec,sf)
         ENDDO
      ENDDO
C
C---        Make sure that the total probability =1 for the whole spectrum in IDSF+1
C
      sfsum=HSUM(IDSF+1)
C      WRITE(6,*) ' A=',nuclaf,' IDSF=',IDSF,' normalization=',sfsum 
      IF(sfsum.GT.0.) CALL HOPERA(IDSF+1,'+',IDSF+1,IDSF+1,1./sfsum,0.) ! normalization
C
C---     IDSF+2 - a smaller area histogram with a finer binning
C
      DO ip=1,npi
         DO ie=1,nei
            dav=0.D0
            DO iavp=1,nav
               DO iave=1,nav
                  ppc=p1+(ip-1.+(iavp-0.5)/REAL(nav))*pbi
                  eec=e1+(ie-1.+(iave-0.5)/REAL(nav))*ebi
                  sf=SF_LOOKUP(ppc,eec)
                  dav=dav+DBLE(sf)
               ENDDO
            ENDDO
            sf=REAL(dav)/nav**2*pbi*ebi
            CALL HF2(IDSF+2,ppc,eec,sf)
         ENDDO
      ENDDO
      IF(sfsum.GT.0.) CALL HOPERA(IDSF+2,'+',IDSF+2,IDSF+2,1./sfsum,0.) ! normalization to the full distribution      
C
      CALL HCOPY(IDSF+1,IDSF+11,'    ')  ! histogram for HRNDM
      CALL HCOPY(IDSF+2,IDSF+12,'    ')  ! histogram for HRNDM
C
      IERR=0
      GO TO 999
C
 997  WRITE(6,*) ' SF_READ: error reading file, line=',nltot
      GO TO 999
C
 998  WRITE(6,*) ' SF_READ: unexpected EOF, file, line=',nltot
      GO TO 999
C
 999  CONTINUE
C      write(6,*) ' error code=',IERR
      END
C
