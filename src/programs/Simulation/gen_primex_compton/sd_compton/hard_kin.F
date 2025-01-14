	subroutine hard_kin
C
	implicit none
C
	include 'kinem.cmn'
	include 'constants.inc'
C
	real*4   rndm
	external rndm
	real*8 theta
	real*8 thetar
	real*8 rgamma
	real*8 Pprime
	real*8 Pprimex
	real*8 Pprimey
	real*8 Pprimez
	real*8 Pgrad
	real*8 Pgradx
	real*8 Pgrady
	real*8 Pgradz
	real*8 Pe
	real*8 Pex
	real*8 Pey
	real*8 Pez
	real*8 thetae
	real*8 relec
	real*8 rgrad
C
	save
C
	theta  = dacos(costheta)
	phi = 2 * pi * rndm(1)
C
	thetar  = dacos(cthgrad)
C
	rgamma = tgtz * dtan(theta)
	xgamma = rgamma * dcos(phi)
	ygamma = rgamma * dsin(phi)
C
	Pprimex = Eprime * dsin(theta) * dcos(phi)
	Pprimey = Eprime * dsin(theta) * dsin(phi)
	Pprimez = Eprime * costheta
	Pprime = dsqrt(Pprimex**2 + Pprimey**2 + Pprimez**2)
C	print * , 'Pprime = ',Pprime
C
	Pgradx = Egrad * dsin(thetar) *  dcos(phi + phirel)
	Pgrady = Egrad * dsin(thetar) *  dsin(phi + phirel)
	Pgradz = Egrad * cthgrad
	Pgrad = dsqrt(Pgradx**2 + Pgrady**2 + Pgradz**2)
C	print * , 'Pgrad  = ',Pgrad
C
C       Account for radiating a hard photon
	Telec = Egamma - Eprime - Egrad
	Eelec = Telec + me
C
	Pe = dsqrt(Eelec * Eelec - me*me)
C
C       Account for radiating a hard photon
	Pez = Egamma - Pprimez - Pgradz
C
	costhe = Pez / Pe
	thetae = dacos(costhe)
C
	Pex =  - Pprimex - Pgradx
	Pey =  - Pprimey - Pgrady
	 
C	print * , 'Pe     = ',dsqrt(Pex**2 + Pey**2 + Pez**2)/Pe
C
C	phie = datan2(Pey,Pex) + pi
	phie = dmod(datan2(Pey,Pex) + 2 * pi, 2 * pi)
C	phie = dmod(phi + pi, 2 * pi)
C
	relec = tgtz * dtan(thetae)
	xelec = relec * dcos(phie)
	yelec = relec * dsin(phie)
C
	dx = xgamma - xelec
	dy = ygamma - yelec
C         
	rmin = dsqrt(dx*dx + dy*dy)
C
	phir = dmod(phi + phirel, 2 * pi)
C
	rgrad = tgtz * dtan(thetar)
	xgrad = rgrad * dcos(phir)
	ygrad = rgrad * dsin(phir)
C
	Pxtot = Pprimex + Pgradx + Pex
	Pytot = Pprimey + Pgrady + Pey
	Pztot = Pprimez + Pgradz + Pez
C
	Pg1(1) = Pprimex
	Pg1(2) = Pprimey
	Pg1(3) = Pprimez
	Pg1(4) = Pprime
C
	Pg2(1) = Pgradx
	Pg2(2) = Pgrady
	Pg2(3) = Pgradz
	Pg2(4) = Pgrad
C
	Pe1(1) = Pex
	Pe1(2) = Pey
	Pe1(3) = Pez
	Pe1(4) = Pe
C
	RETURN
	END
