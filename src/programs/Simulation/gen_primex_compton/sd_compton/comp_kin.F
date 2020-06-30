	subroutine comp_kin
C
	implicit none
C
	include 'kinem.cmn'
	include 'constants.inc'
C
	real*4   rndm
	external rndm
	real*8 theta
	real*8 rgamma
	real*8 Pprime
	real*8 Pprimex
	real*8 Pprimey
	real*8 Pprimez
	real*8 Pe
	real*8 Pex
	real*8 Pey
	real*8 Pez
	real*8 thetae
	real*8 relec
C
	cutg = 0
	cute = 0
C
	theta  = acos(costheta)
	rgamma = tgtz * tan(theta)
	phi = 2 * pi * rndm(1)
	xgamma = rgamma * cos(phi)
	ygamma = rgamma * sin(phi)
C
	Pprimex = Eprime * dsin(theta) * dcos(phi)
	Pprimey = Eprime * dsin(theta) * dsin(phi)
	Pprimez = Eprime * costheta
C
	Telec = Egamma - Eprime
	Eelec = Telec + me
C
	Pe = sqrt(Eelec * Eelec - me*me)
	Pez = Egamma - Eprime*costheta
	Pex = - Pprimex
	Pey = - Pprimey
C
	costhe = Pez / Pe;
	thetae = acos(costhe)
	phie = phi + pi;
	relec = tgtz * tan(thetae)
	xelec = relec * cos(phie)
	yelec = relec * sin(phie)
C
	dx = xgamma - xelec;
	dy = ygamma - yelec;
C         
	rmin = sqrt(dx*dx + dy*dy);
C
	Pg1(1) = Pprimex
	Pg1(2) = Pprimey
	Pg1(3) = Pprimez
	Pg1(4) = Eprime
C
	Pe1(1) = Pex
	Pe1(2) = Pey
	Pe1(3) = Pez
	Pe1(4) = Pe
C
	RETURN
	END
