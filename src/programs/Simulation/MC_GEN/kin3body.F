c*****************************************************************************************
c	Monte Carlo three-body phase space
c
c	R.A.Schumacher, CMU 11-92
c	R.D.Dickson,    CMU 03-2005 - replaces RS version: this one works.
c
c       05-15-15 - make the 4-vectors to have length 5 to add total energy
c
c       This algorithm could be improved to remove the hit-or-miss
c       search for a legal combination of rm23.  Compute the "local limits"
c       given a chosen rm12 and pick uniformly between those limits.
c
	subroutine kin3body(p0,p1,p2,p3,ierr)
	real p0(5),p1(5),p2(5),p3(5),p23(5)
	real p1tmp(5),p2tmp(5),p3tmp(5),p23tmp(5)
	real lmax23,lmin23
	common/monte/iseed
c
	do i=1,3
	p1(i)    = 0
	p2(i)    = 0
	p3(i)    = 0
	p23(i)   = 0
	p1tmp(i) = 0
	p2tmp(i) = 0
	p3tmp(i) = 0
	p23tmp(i)= 0
	end do
c
c	Find absolute kinematic limits of Dalitz variables.
c
	rm0 = p0(4)
	rm1 = p1(4)
	rm2 = p2(4)
	rm3 = p3(4)
c	write(6,*)'>Masses:',rm0,rm1,rm2,rm3
    	amax23 = rm0 - rm1
	amin23 = rm2 + rm3
	amax12 = rm0 - rm3
	amin12 = rm1 + rm2
c	write(6,*)'amax23 > ',amax23
c	write(6,*)'amin23 > ',amin23
c	write(6,*)'amax12 > ',amax12
c	write(6,*)'amin12 > ',amin12
c
c	Pick invariant masses of pairs 1 & 2 and 2 & 3 uniformly in the 
c	allowed range.  This satisfies conditions for "flat" phase space.
c
 10	continue
	rm12 = sqrt(amin12**2 + ((amax12**2 - amin12**2))*(ran(iseed)))
	rm23 = sqrt(amin23**2 + ((amax23**2 - amin23**2))*(ran(iseed)))
c	write(6,*)'rm12 > ',rm12
c	write(6,*)'rm23 > ',rm23
c
c       Find local limits
c
	E2 = (rm12**2 + rm2**2 - rm1**2)/(2*rm12)
	E3 = (rm0**2 - rm12**2 - rm3**2)/(2*rm12)
c	write(6,*)'E2 > ',E2
c	write(6,*)'E3 > ',E3
	lmax23 = sqrt((E2+E3)**2-(sqrt(E2**2-rm2**2)-sqrt(E3**2-rm3**2))**2)
	lmin23 = sqrt((E2+E3)**2-(sqrt(E2**2-rm2**2)+sqrt(E3**2-rm3**2))**2)
c	write(6,*)'lmax23 > ',lmax23
c	write(6,*)'lmin23 > ',lmin23
c
c       Test to see if our point is in the limits
c
	if(rm23.gt.lmax23.or.rm23.lt.lmin23)then
c	   write(6,*)'Trying again...'
	   goto 10
	endif   
c       Point is good, shove it into vectors
	p23(4) = rm23
c
c	Let the decays occur
c
	call kin1r(p0,p1,p23,cos1,pcm,ierr)
	call kin1r(p23,p2,p3,cos2,pcm,ierr)
c	write(6,*)'P0 >',p0
c	write(6,*)'P1 >',p1
c	write(6,*)'P2 >',p2
c	write(6,*)'P3 >',p3
c	write(6,*)'P23>',p23

c	Check momentum conservation

	sumx = p1(1) + p2(1) + p3(1)
	sumy = p1(2) + p2(2) + p3(2)
	sumz = p1(3) + p2(3) + p3(3)
c	write(6,*)'SUM>',sumx,sumy,sumz
	return
	end
