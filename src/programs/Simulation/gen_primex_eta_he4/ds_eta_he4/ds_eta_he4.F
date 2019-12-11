*
* Dec-12-2019, modified by I.Jaegle to accept command line arguments
*
      program ds_eta_he4
      implicit none
*
      double precision PI, TWOPI, PIBY2, DEGRAD
*
      PARAMETER (PI     = 3.141592653589793238d0)
      PARAMETER (TWOPI  = 2.0d0*PI)
      PARAMETER (PIBY2  = 0.5d0*PI)
      PARAMETER (DEGRAD = PI/180.d0)
*
      integer ntetdv, nebin
      parameter (ntetdv = 100, nebin = 1)
*
      character(100) :: num1char
      character(100) :: num2char
      character(100) :: num3char
      character(100) :: num4char
*
      double precision :: etc
      character*(32) :: fnam_coul
      character*(32) :: fnam_stro
      character*(32) :: fnam
*
      real    fc_a1, fc_a2, fs_b1, fs_b2, fs_c1, fs_c2
      common /ff_com/ fc_a1(nebin,ntetdv), fc_a2(nebin,ntetdv),
     &                fs_b1(nebin,ntetdv), fs_b2(nebin,ntetdv),
     &                fs_c1(nebin,ntetdv), fs_c2(nebin,ntetdv)
*
      double precision  ebeam
      real  BCOULM, BNCOHE, BFIINT, BNINCO
      COMMON/pi0_gen_par_com/ ebeam, BCOULM, BNCOHE, BFIINT, BNINCO
*
      integer nterms
      parameter(nterms = 5)
*
      real      thmax
      parameter(thmax = 10.0d0*DEGRAD)
*
* local variables
*
      integer i, j
      double precision dth, th, dom, dsdo(nterms), ds(nterms-1)
*
      integer ic,  jc
      real    ffd1, ffd2, ffd3, ffd4
*
      double precision RECMAS, PI0MAS, U
      parameter(RECMAS = 3.7274177D0) ! He4
      PARAMETER(PI0MAS = 0.5473D0, U = 0.931494D0)
*
      call get_command_argument(1,num1char)
      call get_command_argument(2,num2char)
      call get_command_argument(3,num3char)
      call get_command_argument(4,num4char)
*
      read(num1char,*)etc
      read(num2char,*)fnam_coul
      read(num3char,*)fnam_stro
      read(num4char,*)fnam
*
      open(20, file=fnam_coul,
     &         form='formatted', status='old')
      do i = 1, nebin
        do j = 1, ntetdv
          read(20,'(2(i4,1x),2(1x,E12.5))') ic, jc, ffd1, ffd2
          if(i.ne.ic.or.j.ne.jc) then
            print*, 'Bad data given for coul ff'
            stop
          endif
*
          fc_a1(i,j) = ffd1
          fc_a2(i,j) = ffd2
*
        enddo
      enddo
      close(20)
*
      open(20, file=fnam_stro,
     &         form='formatted', status='old')
      do i = 1, nebin
        do j = 1, ntetdv
          read(20,'(2(i4,1x),4(1x,E12.5))') ic, jc, ffd1, ffd2,
     &                                              ffd3, ffd4
          if(i.ne.ic.or.j.ne.jc) then
            print*, 'Bad data given for stro ff'
            stop
          endif
*
          fs_b1(i,j) = ffd1
          fs_b2(i,j) = ffd2
          fs_c1(i,j) = ffd3
          fs_c2(i,j) = ffd4
*
        enddo
      enddo
      close(20)
*
*
      BCOULM = 510.d0         ! PDG mean
      BNCOHE =   1.d0         ! Roughly
      BFIINT =   1.d0         ! Roughly, from pi0 analysis
      BNINCO =   0.d0         ! Zero for now
*
      dth    =  dble(thmax)/DBLE(ntetdv)
*
      open(23, file=fnam, form='formatted', status='new')
*
      do j = 1, nebin
        ebeam = etc
        do i = 1, ntetdv
          th  = (DBLE(i)-0.5d0)*dth
          CALL CRSEC6(th,dsdo,j,i)
*
          dom = TWOPI * (DCOS(th-dth*0.5d0) - DCOS(th+dth*0.5d0))
*
* Coulomb:
*
          ds(1) = BCOULM*dsdo(1)*dom
*
* Nuclear Coherent:
*
          ds(3) = BNCOHE*dsdo(3)*dom
*
* Their interference:
*
          ds(2) = SQRT(BCOULM*BNCOHE)*dom*
     &           (COS(BFIINT)*dsdo(2)+SIN(BFIINT)*dsdo(5))
*
* Eta incoherent:
*
          ds(4) = BNINCO*dsdo(4)*dom
*
          write(23,'(2(i4,1x),4(1x,E12.5))') j, i, ds
*
        enddo
      enddo
      close(23)
*
      RETURN
      END
*
*     Calculate parametrized cross-section:
*
      SUBROUTINE CRSEC6(thpi0,sumall,ie,it)
      implicit none
*
      integer ie, it
*
      double precision PI, TWOPI, PIBY2, DEGRAD
*
      PARAMETER (PI     = 3.141592653589793238d0)
      PARAMETER (TWOPI  = 2.0d0*PI)
      PARAMETER (PIBY2  = 0.5d0*PI)
      PARAMETER (DEGRAD = PI/180.d0)
*
      double precision  ebeam
      real  BCOULM, BNCOHE, BFIINT, BNINCO
      COMMON/pi0_gen_par_com/ ebeam, BCOULM, BNCOHE, BFIINT, BNINCO
*
      integer ntetdv, nebin
      parameter(ntetdv = 100, nebin = 1)
*
      real    fc_a1, fc_a2, fs_b1, fs_b2, fs_c1, fs_c2
      common /ff_com/ fc_a1(nebin,ntetdv), fc_a2(nebin,ntetdv),
     &                fs_b1(nebin,ntetdv), fs_b2(nebin,ntetdv),
     &                fs_c1(nebin,ntetdv), fs_c2(nebin,ntetdv)
*
      integer nterms
      parameter(nterms = 5)
*
      double precision sumall(nterms), thpi0
*
      double precision aengel
      double precision dsdopr, dsdonc, dsdint, dsdin2, dsdoic,
     &       paul
*
      double precision C1, A1, A2, A3
      double precision TKINRC, t, t1,
     &  Delta2, PPI0DP, EPI0DP, PI0BET, PI0BEM, DISCRI,
     &  TCOULM, TNCOHE
C
      double precision qMeV, qfm, Z, FC1, FC2, FN1, FN2, ARECOI, RECMAS
C
      double precision PI0MAS, U, ALPHA, rgev2mb, Gpi02g, Coul_coef
      PARAMETER(PI0MAS = 0.5473D0, U = 0.931494D0)
      PARAMETER(ALPHA  = 1.D0/137.0359895D0)
      PARAMETER(rgev2mb = 0.38937966D0)           ! GeV^-2 to millibarn
      PARAMETER(Gpi02g = 1.0D-9)                  ! eV to GeV
      PARAMETER(Coul_coef = 8.D0*rgev2mb*ALPHA*Gpi02g)
      parameter(Z = 2.0D0, ARECOI = 4.0D0, RECMAS = 3.7274177D0)      ! He4
*
      double precision ksi, alphae
      parameter(ksi = 0.25d0, alphae = 1.2d0)
*
      double precision  dcosth, dsinth
*
      dcosth = DCOS(thpi0)
      dsinth = DSIN(thpi0)
*
      C1     = (2.D0*ebeam*RECMAS+PI0MAS*PI0MAS)*0.5d0
      A1     = ebeam*ebeam*dcosth**2-(ebeam+RECMAS)**2
      A2     = 2.D0*ebeam*C1*dcosth
      A3     = C1*C1 - PI0MAS*PI0MAS*(ebeam+RECMAS)**2
C
      DISCRI = A2*A2-4.D0*A1*A3
      if(DISCRI.LE.0.D0) then
        print*, 'CRSEC6 ERR: Discriminator is negative'
        stop
      endif
      PPI0DP = (-A2-DSQRT(DISCRI))/(2.D0*A1)
C
      EPI0DP = DSQRT(PPI0DP*PPI0DP+PI0MAS*PI0MAS)
      TKINRC = ebeam-EPI0DP
C
      t      = PI0MAS*PI0MAS-2.D0*(ebeam*EPI0DP-ebeam*PPI0DP*dcosth)
      t1     = (ebeam*dsinth)**2 ! -t ! 4.d0 * ebeam * PPI0DP * dsin(thpi0*0.5d0)**2
C
      Delta2 = -2.D0*RECMAS*TKINRC      ! the same as t in Mandelstamm
C
      qMeV   = 1.D3*DSQRT(2.D0*RECMAS*TKINRC)
      qfm    = qMeV/197.327053D0        ! in 1/fm
C
      FC1 = fc_a1(ie,it)
      FC2 = fc_a2(ie,it)
      FN1 = fs_b1(ie,it) - ksi * fs_c1(ie,it)
      FN2 = fs_b2(ie,it) - ksi * fs_c2(ie,it)
C
      PI0BET = PPI0DP/EPI0DP
      PI0BEM = PI0BET/PI0MAS
C
      TCOULM = Z * DSQRT(Coul_coef*PI0BEM) * PI0BEM
     &           * ebeam**2 * dsinth / ABS(Delta2)
*
* A.Browman et al, Radiative width of the eta meson
* Phys. Rev. Lett vol32, number 19,  13May 1074
*
      TNCOHE = ARECOI * (ebeam**alphae) * DSQRT(16.D0*1.D-3) * dsinth
cc    TNCOHE = ARECOI * (ebeam**alphae) * DSQRT(1.D2*1.D-3) * dsinth
C
      dsdopr = TCOULM*TCOULM*(FC1*FC1+FC2*FC2)        ! in millibarn/sr
      dsdonc = TNCOHE*TNCOHE*(FN1*FN1+FN2*FN2)        ! in millibarn/sr
      dsdint = 2.D0*TCOULM*TNCOHE*(FC1*FN1+FC2*FN2)
      dsdin2 = 2.D0*TCOULM*TNCOHE*(FC2*FN1-FC1*FN2)
C
      sumall(1) = dsdopr
      sumall(2) = dsdint
      sumall(3) = dsdonc
      sumall(4) = 0.d0
      sumall(5) = dsdin2
*
      return
      end
*
