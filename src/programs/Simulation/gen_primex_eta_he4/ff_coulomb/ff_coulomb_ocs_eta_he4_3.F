*
*   Jan-16-2007, I.Larin implementation of
*   S.Gevorgyan calculations for Coulomb formfactor
*   with different nuclear and charge (ocsillator model) densities
*
*   Dec-12-2019, modified by I.Jaegle to accept command line arguments
*
*   f77 ff_coulomb_ocs_eta_he4_3.f -o ff_coulomb_ocs_eta_he4_3 -L $CERN/2004/lib -lpacklib -lmathlib -lkernlib
*
      block data nucl_data
      implicit none
      include 'nucl.inc'
      data  A_n  /4.d0/
      data  a_0  /1.31d0/     ! [fm]
      data  r_p  /0.8d0/      ! [fm]
*
      data  sigma_by2 /0.88d0/  ! eta+N half total c-s at 10GeV in [fm^2]
      data  alpha     /-0.2d0/  ! Re()/Im() elem. amplitude ratio from pi0 to eta
*
      end
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
      parameter (ntetdv =  100, nebin = 1)
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
      dth     =  dble(thmax)/DBLE(ntetdv)
*
      ro_0a   = 1.d0 / (2.d0 * PI * a_0 * a_0)
      rp_a0   = (r_p/a_0)**2 + 1.d0
      ra_0    = dsqrt(r_p**2 + a_0**2)
*
      OPEN(21, FILE=fnam, form='formatted', status='new')
*
      do iebeam = 1, nebin
        ebeam = etc
        do itheta = 1, ntetdv
*
          theta  = (DBLE(itheta)-0.5d0)*dth
*
          call formfactorc(ebeam,theta,w1,w2)
          write(21,'(2(i4,1x),2(1x,E12.5))') iebeam, itheta, w1, w2
          call flush(21)
        enddo
c       call flush(21)
      enddo
*
      close(21)  
      end
*
      subroutine formfactorc(ebeam,theta,w1,w2)
      implicit none
      double precision x
      common /xc/ x(3)
      include 'nucl.inc'
      external  nucl_data
      double precision Q, Q_L
      common/kine/ Q, Q_L
      integer func
      common /pars/ func
*
      external Ft
*
      double precision mpi0, fmgev
      parameter(mpi0 = 0.1349766d0, fmgev = 0.1973269603d0)
      double precision meta
      parameter(meta = 0.5473d0)
*
      double precision ebeam, theta, DGMLT3, W1, W2, CR, CI, c_0
*
      Q_L   =   0.5d0 * meta**2/ebeam /fmgev
      Q     =   2.D0  * ebeam*DSIN(theta*0.5d0) /fmgev
*
      func  = 1
*
      CR   = 1.d0 + 
c
     &       Q*((Q_L/Q)**2+1.d0)*DGMLT3(Ft,0.d0,20.d0,20,8,x)
     &     * 0.5d0/rp_a0**2 * ra_0
*
      func  = 2
*
      CI   = Q*((Q_L/Q)**2+1.d0)*DGMLT3(Ft,0.d0,20.d0,20,8,x)
     &     * 0.5d0/rp_a0**2 * ra_0
*
      c_0   =  dexp((0.5d0*a_0*Q)**2/A_n)
*
	    W1    = c_0 * CR
	    W2    = c_0 * CI
*
      return
      end
*
      subroutine Ft(m,u,F,x)
      implicit none
*
      integer L, m
      double precision u(*), F(*), x(*)
      double precision DGMLT2
      external Fz
*
      do L = 1, m
        x(2)  = u(L)
        F(L)  = DGMLT2(Fz,0,20.D0,20,8,x)
      enddo
*
      return
      end 
*
      subroutine Fz(m,u,F,x)
      implicit none
*
      integer L, m
      double precision u(*), F(*), x(*)
*
      include 'nucl.inc'
      external  nucl_data
      double precision Q,Q_L
      common/kine/ Q,Q_L
      integer func
      common /pars/func
*
      double precision DBESJ1
      double precision b, koe0, koef, fro_p, fro_m, fro_x2m, fro_x2p,
     &  ex_p, ex_m, cql, sql, x22, b2
*
      x22 = x(2)*x(2)
      cql = DCOS(Q_L*x(2)*a_0)
      sql = DSIN(Q_L*x(2)*a_0)
      fro_x2p = sigma_by2*A_n*ro_0a * (1.d0 - derf(x(2)))
      fro_x2m = sigma_by2*A_n*ro_0a * (1.d0 + derf(x(2)))
*
      do L = 1, m
        x(3) = u(L)
           b = DSQRT( (x(3)*x(3)+x22) / rp_a0)
          b2 = b*b
*
       fro_p = fro_x2p*dexp(-x(3)*x(3))
       fro_m = fro_x2m*dexp(-x(3)*x(3))
        ex_p = DEXP(-fro_p)
        ex_m = DEXP(-fro_m)
*
      if (func.eq.1) then
        koef = cql*DCOS(fro_p*alpha)*ex_p
     &       - sql*DSIN(fro_p*alpha)*ex_p
     &       + cql*DCOS(fro_m*alpha)*ex_m
     &       + sql*DSIN(fro_m*alpha)*ex_m
        koe0 = 2.d0*cql
      else if (func.eq.2) then
        koef = cql*DSIN(fro_p*alpha)*ex_p
     &       + sql*DCOS(fro_p*alpha)*ex_p
     &       + cql*DSIN(fro_m*alpha)*ex_m
     &       - sql*DCOS(fro_m*alpha)*ex_m
        koe0 = 0.d0
      endif
*
       F(L)  = DBESJ1(x(3)*Q*a_0) * x(3)**2 / b / b2 *
     & ( koef * ( derf(b) - 2.d0*b/PI05*dexp(-b2) ) - koe0 )
      enddo
*
      return
      end
*
