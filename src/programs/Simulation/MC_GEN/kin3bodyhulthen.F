c******************************************************************************
c	KIN3BODYHULTHEN
c	Monte Carlo three-body phase space
c	Treat it as a quasi-free two-body problem 
c	Take P1 and P2 as the correlated pair which does a (quasi)
c	two-body decay, and take particle 3 as the recoiling particle
c
c	The approach is to keep the recoiling particle (the nucleus) 
c	on shell.  A quasi-particle recoils against the nucleus, and
c	it is off the mass shell.  The incident particle plus the
c	quasi-particle form a system which decays into the other two 
c	final state particles.  
c	
c	This routine has as input the system P0, which consists of the
c	beam particle and the target nucleus.  The recoil nucleus (P3)
c	is allowed to recoil in the lab frame, and the remaining energy 
c	and momentum are computed as P12, which contains the 
c	quasi-particle and the beam particle.  P12 then decay to P1 and P2.
c
c	R.A.Schumacher, CMU 7-95
c	Modifed:
c		2-27-96	- compute fermi-motion recoil in lab frame,
c			  not in beam-plus-target c.m. frame.
c               4-3-99  - noted that in calculating K- + 12C -->
c                         K+ + 11B + Xi- that a small fraction of
c                         events are generated below threshold.  This is
c                         a bug.  A quick review of the code showed no
c                         obvious problems, so it may be related to 
c                         round-off errors.  This needs to be investigated
c                         if this code is to be used in the threshold region.
c               4-13-12 - reuse this code... add some diagnostics
c               9-04-14 - recover code after disk crash and loss of changes since
c                         April of this year
c               9-05-14 - Introduce the Hulthen deuteron momentum distribution for recoil
c               2-27-15 - Put in an optional hack to force a unique value of recoil momentum
c               8-24-18 - repair code after over three years of dormancy
c              10-12-20 - fix /monte/ common block and declare "external ran"
c
c	P0 - c.m. vector for the whole beam-plus-target system in lab frame
c	P1,P2,P3 - as input these contain only the masses
c		   of the three outgoing particles; going out these are the
c		   vectors for the final state particles in the lab frame. 
c
	subroutine kin3bodyhulthen(p0,p1,p2,p3,ierr)
	external ran
	logical ifirst
	real p0(5),p1(5),p2(5),p3(5),p12(5)
	real p1tmp(5),p2tmp(5),p12tmp(5),beta(3)
	data ifirst /.true./
	data phivaluemax /10000.0/
	data alpha / 0.2316 / !inverse Fermi
	data hbarc / .197337/ !GeV fm
	data pmaxgen/0.40/    !GeV/c
c
	common /monte/ iseed,xflat,xnormal,xbw,ixcount
        common/fermi/prec,qfmass,pfermi,ndynamics
	common/hulthen/phivaluemax,alphahc,betahc
c
	do i=1,3
	p1(i)      = 0
	p2(i)      = 0
	p3(i)      = 0
	p12(i)     = 0
	p1tmp(i)   = 0
	p2tmp(i)   = 0
	p12tmp(i)  = 0
	end do
c
c	Find total kinetic energy available in c.m. of particle 0,
c	which is the particle we are going to let decay.
c
c	write(6,*)'Enter kin3bodyhulthen ', p0(4),p1(4),p2(4),p3(4)
	ttot = p0(4) - p1(4) - p2(4) - p3(4)
	if(ttot.lt.0)then
		write(6,*)'No way, no how; masses don''t add up:',
     1  					p0(4),p1(4),p2(4),p3(4)
		return
	end if
c
c	Pick Fermi-motion of target particles uniformly within a sphere...
c
	if(ifirst.eqv..true.)then
	   ifirst=.false.
	   write(6,*)'Max momentum is = ',pmaxgen,' GeV/c'
	   alphahc = alpha*hbarc
	   betahc  = alphahc * 5.98
	endif
c
	imess = 0
100	pxfg = (2.*ran(iseed)-1.)*pmaxgen
	pyfg = (2.*ran(iseed)-1.)*pmaxgen
	pzfg = (2.*ran(iseed)-1.)*pmaxgen
	prec = sqrt(pxfg*pxfg + pyfg*pyfg + pzfg*pzfg)
c
c       Keep this momentum if the overall recoil momemtum is consistent
c       with the deuteron Hulthen momentum distribution.
c
	phivalue = (1./(prec**2. + alphahc**2.) - 
     1              1./(prec**2. + betahc**2))**2.
	ptest = ran(iseed) * phivaluemax
c	write(6,*)prec,phivalue,ptest
	if(phivalue.gt.phivaluemax)then
	   phivaluemax = phivalue
	   write(6,*)'HULTHEN: New largest value of phivalue',
     1               phivaluemax
	endif
c
c       Change the recoil momentum magnitude to be one unique value.
c       This is just a test used to check momentum resolution.
c
cc	pmaxhack = 0.101  !GeV/c
cc	pxfg = pxfg * pmaxhack/prec
cc	pyfg = pyfg * pmaxhack/prec
cc	pzfg = pzfg * pmaxhack/prec
cc	prec = pmaxhack
cc	ptest= 0.0
c
	if(ptest.gt.phivalue)goto 100   !keeps accepted values below curve
c
c	Carve out the recoiling final state particle (3rd one on the list)
c	from the whole system we have been given.  The recoil momentum
c	components we have computed above are already in the lab frame, 
c	so we don't need to boost this particle to the lab.
c
	rmrecoil = p3(4)
	erecoil  = sqrt(prec*prec + rmrecoil*rmrecoil)
	p3(1)  = pxfg 
	p3(2)  = pyfg 
	p3(3)  = pzfg 
	p3(5)  = sqrt(p3(1)**2. + p3(2)**2. + p3(3)**2. + p3(4)**2.)
	p12(1) = p0(1) - pxfg	!Construct the left-over system
	p12(2) = p0(2) - pyfg
	p12(3) = p0(3) - pzfg
	e12    = etot(p0) - erecoil
	rad    = e12**2. - (p12(1)**2. + p12(2)**2. + p12(3)**2.)
	if(rad.gt.0.0)then
		qfmass = sqrt(rad)
	else
		qfmass = 0
	end if
	p12(4) = qfmass
	p12(5) = e12
c
c	Does the kinematic combination make sense?
c
	if(qfmass.lt.(p1(4)+p2(4)))then
	   imess = imess + 1
c	   write(6,*)
c     1	   'Off shell-ness does not work', p12(4)
c	   write(6,*)p0
c	   write(6,*)p1
c	   write(6,*)p2
c	   write(6,*)p3
c	   write(6,*)p12
	   if(imess.gt.1000)then
c	      write(6,*)
c     1         'Off shell-ness is consistently too big', p12(4)
c	      write(6,*)p0
c	      write(6,*)p1
c	      write(6,*)p2
c	      write(6,*)p3
c	      write(6,*)p12
	      ierr = 1
	      return
c	      call exit()
	   else
	      goto 100
	   end if
	end if
c
c	Compute the velocity of the decaying system (quasi-particle plus
c	beam particle) in the lab frame
c
c	betasq = 0
c	e12 = etot(p12)
c	do i=1,3
c	beta(i) = -p12(i)/e12
c	betasq = betasq + beta(i)*beta(i)
c	end do
c
c	Let the decay occur 
c
c	p1tmp(4) = p1(4)
c	p2tmp(4) = p2(4)
c	write(6,*)'Calling     kin1r',p12,p1,p2,cos2
 10	call kin1r(p12,p1,p2,costhcm,pcm,ierr)
c	if(costhcm.gt.0)goto 10  !major hack, for a test

c	write(6,*)'Return from kin1r',p12
c	write(6,*)'                 ',p1
c	write(6,*)'                 ',p2
c	write(6,*)'                 ',costhcm
	if(ierr.ne.0)then
		write(6,*)'KIN1R complains'
	end if
c
c	Check over-all momentum conservation
c
	dpx = p0(1) - p1(1) - p2(1) - p3(1)
	dpy = p0(2) - p1(2) - p2(2) - p3(2)
	dpz = p0(3) - p1(3) - p2(3) - p3(3)
	dp = sqrt(dpx*dpx + dpy*dpy + dpz*dpz)
	if(dp.ge.0.001)then
		write(6,*)'Momentum mismatch:'
		write(6,*)'P0 >',p0
		write(6,*)'P1 >',p1
		write(6,*)'P2 >',p2
		write(6,*)'P3 >',p3
		write(6,*)'P12>',p12,prec
		write(6,*)'dP >',dpx,dpy,dpz,dp
	end if
	return
	end
