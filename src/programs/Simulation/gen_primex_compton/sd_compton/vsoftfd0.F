	real*8 FUNCTION VSOFTFD(X)

	implicit none
	include 'kinem.cmn'
	include 'constants.inc'
C
	real*4   rndm
	external rndm

	real*8 X(2)
	real*8 BORNFD
	real*8 BORNTMP
	external BORNFD
	external BORNTMP
	real H2, P2
	real*8 omega0	! pass energy of photon and cut on the enrgy of emitted photon
C
	real US,UV
	real asinh
C       
	real omega1,y 
	real chi1,chi2
	real ytan
	real cth
C
	real*8 theta
	real*8 rgamma 
C 
        logical gam2in 
C
C       if(costheta.gt.0.9999) costheta=0.9999999

	omega0=X(1)
	costheta=X(2)
	
	omega1 = omega0/(1.+omega0/ME*(1-costheta))
	
	chi1=-2.* omega0/ME
	chi2= 2.* omega1/ME
	
	y=asinh(SQRT(abs(-(chi1+chi2)/4.))) 
	
c       
c       we need to regularize ratio y/atanh(2*y) at y=0.
c       
	if(y.ne.0) then
	   ytan =y/tanh(2*y)
	else
	   ytan = 0.5
	endif
	
	cth=costheta
	
*       soft photon part
	US= BORNTMP(X)*
	1    (4*ytan*(H2(2*y)-1) + 2*(1-2*ytan)*(dlog(2*omega2max/me)-0.5) )
	
*       virtual part
	UV=P2(chi1,chi2,omega1,y,cth,omega0)
	1    +P2(chi2,chi1,omega1,y,cth,omega0)	
	
*       write(*,*)UV,US,X(1),X(2),born(X),H(2*y),y
C
	Egamma = omega0
	Eprime = omega1	
C
	if (.TRUE.) then
	   theta  = dacos(costheta)
	   rgamma = tgtz * dtan(theta)
	   phi = 2 * pi * rndm(1)
	   xgamma = rgamma * dcos(phi)
	   ygamma = rgamma * dsin(phi) 
C      
	   if (.TRUE.) then
	      vsoftcorr = BORNFD(X) - alpha/PI * (US+UV)/omega0
	   else
	      vsoftcorr = 0.d0
	   endif
C       
	endif
C
	VSOFTFD = vsoftcorr
	
	RETURN
	END
