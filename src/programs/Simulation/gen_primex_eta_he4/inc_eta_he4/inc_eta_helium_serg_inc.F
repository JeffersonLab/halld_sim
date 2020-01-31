***************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2007-2019  PrimeEX, GlueX, and PrimEX-D Collaborations     * 
*                                                                         *                                                                                                                               
* Author: The PrimeEX, GlueX, and PrimEX-D Collaboration                  *
* Contributors: Ilya Ilarin, Igal Jaegle                                  *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
***************************************************************************
*
* f77 inc_eta_helium.f -o inc_eta_helium -L $CERN/2004/lib -lpacklib -lmathlib -lkernlib
*
      program main ! gauss,eta,helium 6/03/08
      implicit none
*
      character(100) :: num1char
      character(100) :: num2char
*
      double precision :: etc
      character*(32) :: fnam
*
      double precision x
      common /xc/ x(3)
*
      double precision a_0,  ro_0, A_n, sigma_by2, alpha
      common /nucl/    a_0,  ro_0, A_n, sigma_by2, alpha
*
      double precision Q,Q_L
      common/kine/     Q,Q_L
*
      integer func
      common /pars/ func
*
      external F3
*
      double precision PI, RADDEG
      parameter(PI = 3.141592653589793238d0)
      parameter(RADDEG = PI / 180.d0)
*
      double precision meta, fmgev
      parameter(meta = 0.5475d0, fmgev = 0.1973269603d0)
*
      double precision ebeam, a, b,ap,cp,theta, DGMLT3, W
      real    SR, SI, SS,FF,Wo
      integer n
*
      call get_command_argument(1,num1char)
      call get_command_argument(2,num2char)
*
      read(num1char,*)etc
      read(num2char,*)fnam
*
      ebeam = etc
      Q_L   = 0.5d0*meta**2/ebeam/fmgev
       
       A_n   =  4.d0
      a_0   = 1.292d0   
      ro_0=1.d0/(pi**1.5*a_0**3)

      sigma_by2 =   0.88d0    ! eta+N half total c-s at 10GeV in [fm^2]
      alpha     =  -0.2d0     ! Re()/Im() pi,N amplitude ratio
*
      a     =  0.d0 
      b     = 50.d0
*
      ap    = 0.32d0
      cp    = 2.0d0
*
      open(23, file=fnam, form='formatted', status='new')
*
      do n = 1, 10000
*
        theta = (n-0.1)*0.001
        Q= 2.D0*5.0684*ebeam*DSIN(theta*RADDEG*0.5d0)
*
        func  = 1
        SR    = DGMLT3(F3,a,b,50,8,x)
*
        func  = 2
        SI    = DGMLT3(F3,a,b,50,8,x)
*
        func  = 3
        SS    = DGMLT3(F3,a,b,50,8,x)  
        
        FF    = (2.d0*PI/Q)**2*(SR**2+SI**2)
*
        W=(SS*Q**2 + 0.052d0) * dexp(-ap*Q**2) - A_n*Q**2 * FF
        W= cp * W
*
        Wo    = cp * Q**2 * DEXP(-ap*Q**2) * SS *
     &        ( 1.d0 - (1.d0+0.35d0*Q**4) * DEXP(-1.18d0*Q**2) )
*
        write(23,'(4(E12.6,2x))') theta, W, Wo, SS
      enddo
*
      close(23)
*
      end
*
* ---------------
* Int x(3)=b
*
      subroutine F3(m,u,F,x)
      implicit none
*
      integer m
      double precision u(*), F(*), x(*)
*
      double precision DGMLT2
      external  F1
      integer L
*
      integer func
      common /pars/ func
*
      if(func.eq.1.or.func.eq.2) then
        do L = 1, m
          x(3) = u(L)
          F(L) = DGMLT2(F1,-10.0d0,10.0d0,20,8,x)
        enddo
      else if(func.eq.3) then
        do L = 1, m
          x(3) = u(L)
          F(L) = DGMLT2(F1,  0.0d0, 1.d0,20,8,x)
        enddo
      else
        print*, 'error, wrong func ', func
        call exit(1)
      endif
*
      return
      end 
* ---------------
*Int x(2)=beta
*
      subroutine F1(m,u,F,x)
      implicit none
*
      integer m
      double precision  F(*),u(*),x(*)
*
      double precision a_0, ro_0, A_n, sigma_by2, alpha
      common /nucl/    a_0, ro_0, A_n, sigma_by2, alpha
*
      double precision Q, Q_L
      common/kine/     Q, Q_L
*
      integer func
      common /pars/ func
*
      external Fro
*
      integer L
      double precision  DBESJ1, Fro, frox2, fro15, bes3, dcs,
     & CI, SI, r, pp, abs, sd,ap, cp
*
      
      parameter(ap = 0.32d0, cp = 2.d0)
*
      double precision PI, RADDEG
      parameter(PI = 3.141592653589793238d0)
      parameter(RADDEG = PI / 180.d0)
*
      if(func.eq.1.or.func.eq.2) then
*
        bes3 = DBESJ1(Q*x(3))*x(3)*x(3)
        
        do L = 1, m
          x(2) = u(L)
*
          r  = (x(3)*x(3)+x(2)*x(2))/a_0**2
          PP = 2*ro_0*dexp(-r)/a_0**2
         		
		frox2= sigma_by2*A_n*Fro(x(2))
*
          if(func.eq.1) then
            dcs = DCOS(Q_L*x(2)+alpha*frox2)
          else
            dcs = DSIN(Q_L*x(2)+alpha*frox2)
          endif

          F(L) = (1.d0-sigma_by2*Fro(x(2)))**3*dcs*bes3*PP

        enddo
*
      else if(func.eq.3) then
*
        fro15 = Fro(-15.d0)
*
        do L = 1, m
          x(2) = u(L)
       sd   = 2.d0*sigma_by2 
       abs  = (1.d0 - DEXP(-sd*A_n*fro15) ) / sd
*
       F(L) = 2.d0*pi*x(3)*abs
        enddo
      endif
*
      return     
      end
*
* ---------------
*
      DOUBLE PRECISION FUNCTION Fro(y)
      implicit none
*
      double precision y
      EXTERNAL FS
*
      double precision DGAUSS
*
      Fro = DGAUSS(FS,y,15.d0,1.d-3)
*
      return
      end
*
      DOUBLE PRECISION FUNCTION FS(xg)
      implicit none
*
      double precision  xg
*
      double precision x
      common /xc/ x(3)
*
      double precision a_0, ro_0, A_n, sigma_by2, alpha
      common /nucl/    a_0, ro_0, A_n, sigma_by2, alpha
*
      double precision  r, c43
      parameter(c43 = 1.333333333333333333d0)
*
      r   = (xg*xg + x(3)*x(3))/ a_0**2

      FS = ro_0*dexp(-r)
      return
      end
