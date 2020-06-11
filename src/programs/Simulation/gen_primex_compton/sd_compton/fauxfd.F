
      REAL*8 function BORNTMP(X)
      IMPLICIT NONE
      real*8 X(2)
      REAL*8 costheta
      REAL*8 omega,omega1
      REAL*8 U0,sintheta2
      REAL*8 PI,alpha,me
      parameter (PI = 3.14159 26535 89793 23846 D+00)
      parameter (ALPHA =.00729927 d0)
      parameter (ME = 0.51099892D-3) 
*
      
      omega = X(1)
      costheta = X(2)

      omega1=omega / (1.+ omega / ME * (1 -  costheta))
      sintheta2 = 1 - costheta**2
      U0 = (omega1 / omega + omega / omega1 - sintheta2)

      BORNTMP = PI * ALPHA**2 / ME**2 * (omega1 / omega)**2 * U0
      RETURN
      END

*-----------------------------------------------------------
      real FUNCTION H2 (y)
c     function h(y) from eq.6 PN37
      real y,ABS
      intrinsic abs
      real RDILOG
      real y1
      external RDilog

      real*8 PI,alpha,me
      parameter (PI = 3.14159 26535 89793 23846 D+00)
      parameter (ALPHA =.00729927 d0)
      parameter (ME = 0.51099892D-3) 
      

      if(abs(y) < 0.001) then
         H2 = 1.
         RETURN
      endif
         
      H2=(-Pi**2 - 6 * y**2 + 12 * y * Log(ABS(1. - Exp( 2 * y))) + 
     &     6 * RDIlog(exp(2 * y))) / (12.* y)

		
      RETURN
      END

*------------------------------------------------------------------------

      real FUNCTION P2(hi1,hi2,omega1,y,cth,omega0)
      implicit  none
      real*8 XX(2)
      real y, log,abs,tanh,cosh,sinh,hi1,hi2,Omegasv1,Omega1,RdILOG
      real H2, u0sv
      External U0sv,RdILOG
      intrinsic  log,abs,tanh,cosh,sinh
      real cth
      real ytan2, ytan
      real U0
      real as

      real*8 X
      real*8 omega0
      real*8 BORNTMP,costheta
      REAL*8 PI,alpha,me
      parameter (PI = 3.14159 26535 89793 23846 D+00)
      parameter (ALPHA =.00729927 d0)
      parameter (ME = 0.51099892D-3) 

      if(y.NE.0) then
         ytan2 = y * (1 / tanh(2 * y)) 
         ytan =y * (1 / tanh(y)) 
      else
         ytan2 = 0.5
         ytan = 1
      endif

      if(abs(hi1 + hi2).GT.0.0001) then
         as = y**2 / (hi1 + hi2)
      else 
         as=-0.25
      endif


      costheta = cth
      XX(1) = omega0
      XX(2) = costheta


      U0 = BORNTMP(XX) / pi / alpha**2 * ME**2 * omega0**2 / omega1**2


      P2 = -12/hi1 + 4*(1/hi1 + 1/hi2)**2 + (0.5+hi1/hi2)/(hi1-1) -
     &     (2*hi1)/hi2**2 - (3*hi1)/(2.*hi2)- 
     -     4*(2/hi1-(7*hi1)/4.-(3*hi2**2)/(4.*hi1))*as+ 
     -     (2*(-1 + 2/hi1 + hi1 - 3/hi2 + hi1/hi2 + hi1**2/hi2 + hi2/2.+ 
     &     hi2/hi1**2)*(RDIlog(1.) - RDIlog(1. - hi1)))/hi1
     -     - 2*U0*ytan2*(2*H2(y) - H2(2*y)) + 
     -     Log(Abs(hi1))*(1 - 8/hi1**2 + 8/hi1 + 3/hi2 - 7/(hi1*hi2) + 
     &     (3*hi2)/(2.*hi1**2) + (3*hi2)/(2.*hi1) - 
     -     (2*hi1**2 + hi2)/(2.*(-1 + hi1)**2*hi2) + (2*hi1 - hi1**2*hi2-
     &     hi2**2)/(2.*(-1 + hi1)*hi1**2*hi2) + 
     -     4*ytan2*(-1 + 4/hi1**2 - 1/hi1 - hi1/hi2 - 
     &     hi2/(2.*hi1) + (4*Cosh(y)**2)/(hi1*hi2) + 
     -     ((-6 + hi1)*(1/Cosh(2*y)))/(2.*hi2))) 
     &     + H2(y)*(2*ytan-(4*y*(2 - Cosh(2*y))*Sinh(2*y))/(hi1*hi2)) 
     &     -4*(0.5 - 1/hi1)*y*Tanh(y)
	

      P2 = PI * ALPHA**2 / ME**2 * (omega1 / omega0)**2 * P2
	
      RETURN
      END
