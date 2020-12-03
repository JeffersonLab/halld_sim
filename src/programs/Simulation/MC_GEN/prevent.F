c*****************************************************************************
c	PREVENT
c
c	Prints out the details of one event block
c
c	Coding for the list of particles generated in the decay chain:
c	All variables in this list are measured in the lab.
c	PLIST(I,MARK)	
c	I = index into list of a given particle's parameters
c	MARK = counts the particles produced
c	I = 
c	1	PX
c	2	PY
c	3	PZ
c	4	MASS
c	5	TOTAL ENERGY
c	6	LUND ID FOR THIS PARTICLE
c	7	X OF CREATION
c	8	Y OF CREATION
c	9	Z OF CREATION
c	10	X COMPONENT OF POLARIZATION
c	11	Y COMPONENT OF POLARIZATION
c	12	Z COMPONENT OF POLARIZATION
c	13	PX AT DECAY POINT
c	14	PY AT DECAY POINT
c	15	PZ AT DECAY POINT
c	16	MASS
c	17	TOTAL ENERGY AT DECAY POINT
c	18	X OF DECAY
c	19	Y OF DECAY
c	20	Z OF DECAY
c	21	CTAU
c
c	This subroutine made into its own file 7-17-01
c       9-23-2014 - tweak output to give time-of-flight over 4 meters.
c       3-23-2016 - add beam momentum to output
c
	subroutine prevent(plist,inext,ievent)
	parameter (maxnumlis=50) !Maximum number of particles in the decay list
	real plist(24,maxnumlis)
	common /beam/pbeam,ppbeam,rmbeam,rmtarg,jcode,xsigma,ysigma,ncode
	write(6,*)'Event, Beam momentum',ievent,ppbeam
	do i=1,inext-1
	   beta = absp(plist(1,i))/plist(5,i)   !pc/E
	   tof  = 4.0/(beta*0.3) !nanoseconds over 4 meter flight path
	   write(6,10)int(plist(6,i)),(plist(j,i),j=1,5),  beta,tof,
     1     (plist(j,i),j=13,15),
     1     plist(17,i),
     1     (plist(j,i),j=7,9),(plist(j,i),j=10,12),(plist(j,i),j=18,20)
 10	   format(1x,i4,') px,py,pz,m,E ',5f10.3,
     1     '            beta, tof ',f10.3,f10.3,/,
     1     6x,'              ',3f10.3,10x,f10.3,/
     1     6x,'X,Y,Z,Polx,y,z',3(1pe10.1),3f10.4,/,
     1     6x,'              ',3(1pe10.1))
c	   call clascrunch(int(plist(6,i)),plist(1,i),plist(2,i),
c     1	                       plist(3,i), plist(4,i),plist(5,i))
	enddo
	return
	end

c
c       CLAS-specific code, hacked together from pieces of
c       mc_analyze.
c
c       Feb. 2007
c
	subroutine clascrunch(lund,px,py,pz,pm,Etot)
	data rade /57.29577951/
c
	dummy = pm  !to fool compiler
	dummy = Etot
	pmom = sqrt(px*px +py*py +pz*pz) !MeV/c
	if(pmom.gt.0)then
		arg=pz/pmom
		if(arg.gt. 1.0)arg= 1.0
		if(arg.lt.-1.0)arg=-1.0
		theta_polar = acos(arg) * rade
		phi = atan2(py,px) * rade
		if(phi.lt.0.0)phi = phi + 360.  !azimuthal angle from 0 to 360.
	else
		theta_polar = 1000.
		phi = 1000.
	end if	
	pmomgev = pmom/1000.            !particle momentum in GeV/c 
	phispaz = mod(phi+30.,60.)-30.   !phi within one sector +- 30 deg.
	isector = int((phi+30.)/60.) + 1
	if(isector.eq.7)isector = 1
	if((lund.eq.-17).or.(lund.eq.+7).or.(lund.eq.-18))then
		icharge = -1
	else if((lund.eq.17).or.(lund.eq.18).or.(lund.eq.41))then
		icharge = +1
	else
		icharge = 0
	end if
c	if(icharge.ne.0)then
c	   write(6,10)icharge,isector,theta_polar,phi,phispaz
c 10	 format(1x,'      CLAS: charge, sector, theta, phi, phi_sect ',
c     1         2i5,3f10.1)
c	endif
	return
	end
