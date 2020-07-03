C==========================================================================C
C                                                                          C 
C Hard photon emission contribution ot Compton cross section               C
C Input is a 4-dimentional array XX                                        C
C                                                                          C
C XX(1) = beta1   -- transformation of cos(theta1) angle of first photon   C
C   the costheta1 is defined by using XX(1)=beta1 variable                 C
C XX(2) = beta2   -- transformation of cos(theta1) angle of second photon  C
C  the costheta2 is defined by using XX(2)=beta2 variable                  C
C XX(3) = phi     -- azimuthal angles between two photons                  C
C XX(4) = omega1  -- energy of first photon                                C
C                                                                          C
C The omega0 and omega2max has to be passed by common block                C
C           Common/par/Omega0,omega2max                                    C
C                                                                          C
C Limits of integrations:      a(i)<XX(i)<b(i)                             C
C                                                                          C
C  Here we define a(i)'s and b(i)'s                                        C
C                                                                          C
C--------  for XX(1)                                                       C
C                                                                          C
C      beta=dSqrt(1.d0-2.d0*me/Omega0)                                     C
C      a(1)=dLog((1. - 0.99*beta)/(1. + 0.99*beta))/(2.*beta)              C
C      b(1)=dLog((1. + beta)/(1. - beta))/(2.*beta)                        C
C                                                                          C
C-------- for XX(2)                                                        C
C                                                                          C
C      a(2)=dLog((1. - 0.99*beta)/(1. + 0.99*beta))/(2.*beta)              C
C      b(2)=dLog((1. + beta)/(1. - beta))/(2.*beta)                        C
C                                                                          C
C-------- for XX(3)                                                        C
C                                                                          C
C      a(3)=0.d0                                                           C
C      b(3)=2.*pi                                                          C
C                                                                          C
C------- for XX(4)                                                         C
C                                                                          C
C      beta=dSqrt(1.d0-2.d0*me/Omega0)                                     C
C      dc1=1./beta*dtanh(beta*X(1))                                        C
C      dc2=1./beta*dtanh(beta*X(2))                                        C
C      omegas = omega0/(1.+omega0/me*(1.-dc1))                             C
C      costheta12=dc1*dc2+sqrt(1-dc1*dc1)*sqrt(1-dc2*dc2)*cos(x(3))        C
C                                                                          C
C      f12 = ( (1+omega0/me*(1-dc1))*(1+omega0/me*(1-dc2))                 C
C     & -omega0/me*(1-costheta12))                                         C
C     &/(1+omega0/me*(1-dc1))                                              C
C     &/(1+omega0/me*(1-dc1)-(1-costheta12)*omega2max/me)                  C
C                                                                          C
C      omega1max = omegas - omega2max*f12                                  C
C      a(4)=omega2max                                                      C
C      b(4)=omega1max                                                      C
C                                                                          C
C==========================================================================C
      	
	real*8 function hardfd(XX)
C       
	Implicit none
C       
	include 'kinem.cmn'
	include 'constants.inc'
C
	real*4   rndm
	external rndm
C       
	real*8 XX(5)
	Real*8 a1, a2, a3, b1, b2, b3 ! functions from eq.13 in the ref.4 PN37
	Real*8 dc1,dc2,fi1,fi2
	Real*8 Omega2		!energy of the second photon
	Real*8 e2		!energy of the scattering electron
	
C	Real*8 cosfie,costhe,sinfie1 !cosfie=Cos(fie),costhe=Cos(Thetae),sinfie1=Sin(fie)
	Real*8 L, T1		!functions from eq.15 in the ref.4 PN37
C	Real*8 fie,Thetae	!azimuthal and polar angle of the scattering electron
	Real*8 fi		!fi=fi2-fi1 
	
	
C	Real*8 k1(4),k2(4),p2(4)
C	Real*8 k1x,k1y,k1z,k2x,k2y,k2z,p2x,p2y,p2z ! x,y,z component of the 4-momenta k1,k2,p
C	Real*8 px,py,pz,e	! variables to verify the conservation of momentum and energy in double Compton scattering 
	
	Real*8 Omega0,Omega1,dsigma2
	Real*8 X		!function from eq.9 in the ref.4 PN37
	real*8 costheta12,f12
	real*8 coeff
	real*8 alpha1, beta,gamma,delta,rho,y,A,B
	real*8 beta1
	
	real*8 dc1min, dc1max
	real*8 dc2min, dc2max
	real*8 xx2min, xx2max
	real*8 xx3min, xx3max
	real*8 xx2, xx3
	real*8 phimin, phimax
	real*8 omega1min, omega1max
C
	real*8 theta
	real*8 rgamma
C
	real*8 thgrad
	real*8 rgrad
C
	real*8 Pprimex
	real*8 Pprimey
	real*8 Pprimez
	real*8 Pgradx
	real*8 Pgrady
	real*8 Pgradz
	real*8 Pe
	real*8 Pex
	real*8 Pey
	real*8 Pez
	real*8 thetae
	real*8 relec
C
        logical gam1in, gam2in, eminin
C	
	real*8 omegas
	common/omegas/omegas,f12,costheta12
	
	fi = XX(4)
	Omega0=XX(5)
	
	beta1=dSqrt(1.d0-2.d0*me/Omega0)
	
	xx2min=dLog((1. - 0.99*beta1)/(1. + 0.99*beta1))/(2.*beta1)
	xx2max=dLog((1. + beta1)/(1. - beta1))/(2.*beta1)
	
	xx3min=dLog((1. - 0.99*beta1)/(1. + 0.99*beta1))/(2.*beta1)
	xx3max=dLog((1. + beta1)/(1. - beta1))/(2.*beta1)
	
	xx2=xx2min+(xx2max-xx2min)*XX(2)     
	xx3=xx3min+(xx3max-xx3min)*XX(3) 
	
	dc1=1./beta1*dtanh(beta1*xx2)
	dc2=1./beta1*dtanh(beta1*xx3)
	
	omegas = omega0/(1.+omega0/me*(1.-dc1))
	costheta12=dc1*dc2+sqrt(1-dc1*dc1)*sqrt(1-dc2*dc2)*cos(fi)
	
	f12 = ( (1+omega0/me*(1-dc1))*(1+omega0/me*(1-dc2))
	1    -omega0/me*(1-costheta12))
	2    /(1+omega0/me*(1-dc1))
	3    /(1+omega0/me*(1-dc1)-(1-costheta12)*omega2max/me)
	
	omega1min = omega2max
	omega1max = omegas - omega2max*f12
	omega1 = omega1min + (omega1max - omega1min) * XX(1)
	
	L=me*Omega0 - Omega1*(me + Omega0*(1.D0 - dc1))
	
	T1=me + (1.D0 - dc2)*Omega0 - 
	1    Omega1*(1.D0 - dc1*dc2 - 
	2    dsqrt(1.D0-dc1**2)*dsqrt(1.D0-dc2**2)*dCos(fi))
	Omega2=L/T1
	
	
	a1= Omega1/me
	a2= Omega2/me
	a3= -Omega0/me
	b1= -((-(dc1*Omega0*Omega1) + Omega1*(me + Omega0 - 
	1    Omega2) + Omega1*Omega2*(dc1*dc2 + 
	2    dsqrt(1.-dc1**2)*dsqrt(1.-dc2**2)*dCos(fi)))/me**2)
	b2=(-((-1. + dc1)*Omega0*Omega1) + 
	1    me*(-Omega0 + Omega1))/me**2
	b3=(dc1*Omega0*Omega1 + Omega0*(me - Omega1 - Omega2) + 
	1    dc2*Omega0*Omega2)/me**2
	
	
	
	alpha1 =1./a1+1./a2+1./a3
	beta  =1./b1+1./b2+1./b3
	gamma=1./(a1*b1)+1./(a2*b2)+1./(a3*b3)
	A=a1*a2*a3
	B=b1*b2*b3
	delta=a1*b1 + a2*b2 + a3*b3
	rho=a1/b1+b1/a1  + a2/b2+b2/a2  +  a3/b3+b3/a3
	y=a1+a2+a3
	
	
	X = -2*(alpha1*beta-gamma)**2+4*(alpha1*beta-gamma)*(alpha1+beta)
	1    -16.*alpha1*beta+8.*gamma-2*y*(alpha1**2+beta**2)
	2    +2*(alpha1*beta-gamma)*y*(alpha1+beta)
	3    -4.*y/A/B*(y**2*(delta-1.)-2.*delta)+4*y*(y+1)*(1/A+1/B)
	4    -4*(2.*y+delta*(1.-y))*(alpha1/B+beta/A)
	5    - 2*(alpha1*beta+gamma*(1-y))*rho
	
	
*       
	coeff = 1./(4.*pi)*alpha**3/ME**3/4.
	dsigma2=coeff*(X*Omega1*Omega2/Omega0/T1)/Omega0
	
	if(dsigma2.ne.dsigma2) then
	   write(*,*)'dsigma2',dsigma2
	   stop
	endif
C
C       Fill in KINEM common block
C
	Egamma = Omega0
	Egrad  = Omega1
	Eprime = Omega2
	phirel = fi
	cthrel = costheta12
	cthgrad = dc1
	costheta = dc2	
C
        theta  = dacos(costheta)
        rgamma = tgtz * dtan(theta)
        phi = 2 * pi * rndm(1)
        xgamma = rgamma * dcos(phi)
        ygamma = rgamma * dsin(phi)
C     
        thgrad  = dacos(cthgrad)
        rgrad = tgtz * dtan(thgrad)
        xgrad = rgrad * dcos(phi + phirel)
        ygrad = rgrad * dsin(phi + phirel)
C
	Eelec = Egamma - Eprime - Egrad + me
	Pe = dsqrt(Eelec * Eelec - me*me)
	Pez = Egamma - Pprimez - Pgradz
	costhe = Pez / Pe
	thetae = dacos(costhe)
	Pex =  - Pprimex - Pgradx
	Pey =  - Pprimey - Pgrady
	phie = dmod(datan2(Pey,Pex) + pi, 2 * pi)
	relec = tgtz * dtan(thetae)
	xelec = relec * dcos(phie)
	yelec = relec * dsin(phie)
C     
	if (.TRUE.) then
           hardcorr = dsigma2
     $          *(1-(dtanh(beta1*xx2))**2)
     $          *(1-(dtanh(beta1*xx3))**2)
     $          *(xx2max - xx2min)
     $          *(xx3max - xx3min)
     $          *(omega1max - omega1min)
        else
           hardcorr = 0.d0
        endif
        
        hardfd = hardcorr

	return
	
	END
