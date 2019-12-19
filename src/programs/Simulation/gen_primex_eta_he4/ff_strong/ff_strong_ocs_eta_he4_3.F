***************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2007-2019  PrimeEX, GlueX, and PrimEX-D Collaborations     * 
*                                                                         *                                                                                                                               
* Author: The PrimeEX, GlueX, and PrimEX-D Collaboration                  *
* Contributors: Ilya Ilarin, S. Gevorgyan, Igal Jaegle                    *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
***************************************************************************
*
*   Jan-17-2007, I.Larin implementation of
*   S.Gevorgyan calculations for strong formfactor
*   with finit nucleus radius and ocsillator model density
*
*   Dec-12-2019, modified by I.Jaegle to accept command line arguments
*
*   f77 ff_strong_ocs_eta_he4_3.f -o ff_strong_ocs_eta_he4_3 -L $CERN/2006/lib -lpacklib -lmathlib -lkernlib
*
      block data nucl_data
      implicit none
      include 'nucl.inc'
      data  A_n  /4.d0/
      data  a_0  /1.31d0/     ! [fm]
      data  r_p  /0.8d0/      ! [fm]
*
      data  sigma_by2 /0.88d0/  ! pi0+N half total c-s at 5GeV in [fm^2]
      data  alpha     /-0.2d0/  ! Re()/Im() elem. amplitude ratio
*
      end
*
*
      program main
      implicit none
*
      character(100) :: num1char
      character(100) :: num2char
*
      double precision :: etc
      character*(32) :: fnam
*
      integer ntetdv, nebin
      parameter (ntetdv =  10000, nebin = 1)
*
      include 'nucl.inc'
      external  nucl_data
*
      real      thmax
      parameter(thmax = 10.0*DEGRAD)
*
      integer iebeam, itheta
      double precision ebeam, theta, dth, w1, w2
*
      call get_command_argument(1,num1char)
      call get_command_argument(2,num2char)
*
      read(num1char,*)etc
      read(num2char,*)fnam
*
      dth    =  dble(thmax)/DBLE(ntetdv)
*
      OPEN(21, FILE=fnam, form='formatted', status='new')
*
      do iebeam = 1, nebin
        ebeam = etc
        do itheta = 1, ntetdv
*
          theta  = (DBLE(itheta)-0.5d0)*dth
*
          call formfactors(ebeam,theta,w1,w2)
          write(21,'(2(i6,1x),4(1x,E12.5))')
     &          iebeam, itheta, w1, w2, 0., 0.
cc        call flush(21)
        enddo
        call flush(21)
      enddo
*
      close(21)  
      end
*
      subroutine formfactors(ebeam,theta,w1,w2)
      implicit none
      double precision x(5)
*
      include 'nucl.inc'
      external  nucl_data
*
      double precision Q, Q_L, sl
      common /kine/    Q, Q_L, sl
      integer func
      common /pars/ func
*
      external F3
*
      double precision ebeam, a, b, theta, c_0, DGMLT3,
     &  W1, W2, CR, CI
*
      double precision mpi0, fmgev
      parameter(mpi0 = 0.1349766d0, fmgev = 0.1973269603d0)
      double precision meta
      parameter(meta = 0.5473d0)
*
      sl    = 0.24d0
*
      Q_L   =   0.5d0 * meta**2/ebeam /fmgev
      Q     =   2.D0  * ebeam*DSIN(theta*0.5d0) /fmgev
*
      a     =  0.d0
      b     =  7.d0
*
      func  = 1
      CR    = DGMLT3(F3,a,b,5,8,x)
*
      func  = 2
      CI    = DGMLT3(F3,a,b,5,8,x)
*
	    c_0   = 2.d0 * dexp((0.5d0*a_0*Q)**2/A_n) /sl**2 /Q /PI05 /a_0**3
	    W1    = c_0 * CR
	    W2    = c_0 * CI
*
      return
      end
*
* ---------------
* Int x(3)=b
*
      subroutine F3(m,u,F,x)
      implicit none
*
      integer L, m
      double precision u(*), F(*), x(*)
*
      double precision Q, Q_L, sl
      common /kine/    Q, Q_L, sl
*
      double precision DGMLT2
      external F2
*
      do L = 1, m
        x(3)  = u(L)
        F(L)  = DBESJ1(x(3)*Q)*x(3) * DGMLT2(F2,-5.d0,5.d0,5,8,x)
      enddo
*
      return
      end 
*
* ---------------
* Int x(2)=z
*
      subroutine F2(m,u,F,x)
      implicit none
*
      integer L, m
      double precision u(*), F(*), x(*)
      double precision b, DGMLT1
*
      external F1
*
      integer func
      common /pars/ func
*
      do L = 1, m
        x(2) = u(L)
        F(L) = DGMLT1(F1,0.d0,10.d0,5,8,x)
      enddo
*
      return
      end
*
* ---------------
* Int x(1)=x
*
      subroutine F1(m,u,F,x)
      implicit none
*
      integer L, m
      double precision  u(*), F(*), x(*)
*
      include 'nucl.inc'
      external  nucl_data
*
      double precision Q, Q_L, sl
      common /kine/    Q, Q_L, sl
*
      integer func
      common /pars/ func
*
      double precision DBESJ1, DBESI0, DBESI1, FRO
      double precision a, koef, dcs
*
      double complex   cabsorp, ii
*
      ii   = DCMPLX(0.d0,1.d0)
*
      FRO  = sl/(a_0*a_0+2.d0*sl) *
     &       dexp(-x(3)**2/(a_0*a_0+2.d0*sl)) *
     &       0.5d0/PI * (1.d0 - derf(x(2)/a_0))
*
      cabsorp = (1.d0-sigma_by2*FRO*(1.d0-alpha*ii)/sl)**(A_n-1.d0)
     &        * cdexp(ii*Q_L*x(2))
*
      if (func.eq.1) then
        dcs   = DREAL(cabsorp)
      else if (func.eq.2) then
        dcs   = DIMAG(cabsorp)
      endif
*
      do  L = 1, m
        x(1)  =  u(L)
        koef  = x(3) * DBESI0(x(3)*x(1)/sl) -
     &          x(1) * DBESI1(x(3)*x(1)/sl)
*
        a     = DEXP( -(x(3)*x(3)+x(1)*x(1))/(2.d0*sl)
     &                -(x(1)*x(1)+x(2)*x(2))/a_0**2 )
        F(L)  = x(1) * koef * a * dcs
      enddo
*
      return     
      end
*
