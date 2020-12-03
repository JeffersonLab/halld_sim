
	REAL*8 function BORNFD(X)
	IMPLICIT NONE 
        include 'kinem.cmn' 
        include 'constants.inc' 
C 
        real*4   rndm 
        external rndm 

	real*8 X(2)
	REAL*8 omega,omega1
	REAL*8 U0,sintheta2
C 
        real*8 theta 
        real*8 rgamma  
C  
        logical gam2in  
C
	omega = X(1)
	costheta = X(2)
C
	omega1 = omega / (1.+ omega / ME * (1 - costheta))
	sintheta2 = 1 - costheta**2
	U0 = ( omega1 / omega + omega / omega1 - sintheta2)
C
	if (.TRUE.) then 
           theta  = dacos(costheta) 
           rgamma = tgtz * dtan(theta) 
           phi = 2 * pi * rndm(1) 
           xgamma = rgamma * dcos(phi) 
           ygamma = rgamma * dsin(phi)  
C
           gam2in = dabs(datan(sqrt(xgamma**2+ygamma**2)/tgtz)*180/PI-thctr)<thbin  
C   
C           if (gam2in) then 
           if (.TRUE.) then 

	      BORNFD = PI * ALPHA**2 / ME**2 * (omega1 / omega)**2 * U0/omega
	   else
	      BORNFD = 0.0
	   endif
	endif
C
	RETURN
	END
