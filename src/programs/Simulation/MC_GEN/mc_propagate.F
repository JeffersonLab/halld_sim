c******************************************************************************
c	MC_PROPAGATE
c	Routine to propagate a particle until it decays or stops.
c       This version assumes a uniform target material surrounded by
c       a shell of "wall" material.
c
c       The version of this code for CLAS computed target energy loss.  But this
c       version has had had the energy loss calculation stripped out.  It could probably
c       be further trimmed down.
c
c	R.A.Schumacher, CMU, 3-96
c	Modified:	3-5-96 	- check whether particle leaves target.
c			3-7-96	- better logic for checking target geometry
c                       8-19-97 - UNIX version
c                       5-20-99 - mods for cylindrical (CLAS) target
c                       1-09-14 - test final location r1 for NaN  
c                       5-15-15 - make PV have 5 elements, with energy in element 5
c                       7-05-18 - resurrect this code for GlueX Lambda's;  no energy losses calculated
c                       6-26-19 - slim code down to remove energy loss calculation
c                      10-12-20 - updte /monte/ common block and declare "external ran"  	
c       
c	PV - vector for particle: px, py, pz, mass
c	   - used both for input and output
c	R0 - starting point in X,Y,Z coordinates                -input
c	R1 - calculated decay point in X,Y,Z coordinates        -output
c	CTAU - c-tau value of decaying particle in centimeters  -input
c	ZPAR - charge of particle in units of elementary charge -input
c
	subroutine propagate(pv,r0,r1,ctau,zpar,lund)
	external ran
	real pv(5),r0(3),r1(3),direc(3),r0tmp(3),r1tmp(3)
	real pvtmp(5)
	common /target/icode,size1,size2,size3,size4,maxtstep
	common /dedxlocal/ztarg,zda,riadj,rho
	common /monte/ iseed,xflat,xnormal,xbw,ixcount
	common /temp/eloss_clas,eloss_prop,pcmag_clas,pvmag_prop
c
	pvmag = absp(pv)
	et    = etot(pv)
	rmass = pv(4)
	do i=1,5
	   pvtmp(i) = pv(i)
	enddo
	eloss_prop = 0.0
	pvmag_prop = pvmag   !for later plotting via common
c
c	direction cosines
c
	do kk=1,3
	if(pvmag.gt.0)then
		direc(kk)=pv(kk)/pvmag
	else
		direc(kk)=0.0
	end if
	enddo
	theta = acos(direc(3))   !angle of track in radians with respect to beam axis
c
	r0tmp(1) = r0(1)
	r0tmp(2) = r0(2)
	r0tmp(3) = r0(3)
c
c	For all particles we propagate in one step.
c	This ignores neutron interactions in the target as well as 
c	photon interactions.
c
	tdecay = -ctau*alog(ran(iseed)) !decay "time" in particle rest frame
	deltat = tdecay / maxtstep 	!size of rest frame "time" steps 
	tpart  = et - rmass	!particle kinetic energy in GeV
	tpartstart = tpart      !initial K.E. in GeV
	pvmagstart = pvmag
c	write(6,*)' '
c	write(6,'(1x,a,f10.3,e10.3,f10.3,1pe10.3)')
c     1     'rmass,ctau,tpart,tdecay',rmass,ctau,tpart*1000.,tdecay
	iflag = -1
c
	pvmag = absp(pv)
	et    = etot(pv)
	if(et.le.0)then
		write(6,*)'Bogus particle',pv,r0,ctau,zpart
		write(6,*)' '
		r1(1) = r0(1)
		r1(2) = r0(2)
		r1(3) = r0(3)
		return
	end if
	beta  = pvmag/et
	if(beta.lt.1.0)then
		gamma = sqrt(1./(1.-beta*beta))
	else
		gamma = 100000
	end if
c
c	Take one giant step
c
	if(et.gt.0)then
	   dstep = beta*gamma*tdecay !(remaining) lab distance
	else
	   dstep = 1.e32
	end if
	call travel(r0tmp,direc,dstep,r1)
	eloss_prop =   sqrt(pvtmp(1)**2.+pvtmp(2)**2.+pvtmp(3)**2.)
     1  	-  sqrt(pv(1)**2.   +pv(2)**2.+   pv(3)**2.)
c       write(6,*)'Propagating',lund,zpar,iflag,icode
c       write(6,*)r0
c       write(6,*)r1
	return
	end


c       
c 	FUNCTION DEDX (HEAVY PARTICLES)  CALCULATE DE/DX IN MEV/GCM-2
c 	Calculates dE/dx based on NASA report SP-3013, Barkas and Berger
c
c	RMASS and TPART are in GeV, result DEDX is in MeV
c
c	Modified:
c	3-10-96	- rationalize calls and add more terms
c
	function dedx(zpart,rmass,tpart)
	parameter (rmelec=.511 , twormelec=1.022)	!electron mass (MeV)
	parameter (rk=.153538)				!MeV/gm/cm2
c	real*8 gamma,beta,eta,r,wmax,cadj,arg,delta
	real   gamma,beta,eta,r,wmax,cadj,arg,delta
	common /dedxlocal/ztarg,zda,riadj,rho
c
	dedx = 0.
	if(tpart.eq.0)return
c	if(riadj.eq.0)riadj = 16. * ztarg**0.9		!in eV
 	gamma = (tpart+rmass)/rmass
	beta  = sqrt(1. - 1./gamma**2.)
	eta   = gamma*beta
	r     = (rmelec*0.001)/rmass		!electron mass to particle mass
	wmax  = twormelec * eta**2. / (1. + 2.*r*gamma + r**2.)	!in MeV
	arg   = twormelec * eta**2. * wmax / (riadj*1.e-6)**2.	!unitless
c
c	Shell correction
c
 	cadj  = (.4224/eta**2.+.0304/eta**4-.00038/eta**6.)*riadj**2.*1.e-6
     1       + (3.858/eta**2.-.1668/eta**4+.0016 /eta**6.)*riadj**3.*1.e-9
c
c	Density effect
c
	if(eta.gt.3.0)then
		homega = 28.816*sqrt(rho*zda)		!in eV (PDG '94, p1251)
		delta = 2.*(log(homega/riadj) + log(eta) - 0.5)
	else
		delta = 0.0
	end if
c
	dedx  = rk * zpart**2. * zda /beta**2. * 
     1  		(log(arg) - 2.*beta**2. - 2. * cadj/ztarg - delta)
c	write(6,*)'dedx',zpart,rmass,tpart,dedx
c	write(6,*)gamma,beta,eta,r,wmax,arg,cadj,delta
c	write(6,*)ztarg,zda,riadj
	return
	end

c
c	Routine to see if a given step takes us across a target boundary.
c	If so, reset the destination point to be ON the boundary.
c	The logic is for a cylindrical target aligned with the coordinate
c	axes.  In this version the cyl. is also centered on the origin
c
	subroutine checkout(r0,r1,direc,deltad,iflag)
	logical ifirst
	real r0(3),r1(3),direc(3)
	common /target/icode,size1,size2,size3,size4,maxtstep
	common /checklocal/radius1,radius2,zmin,zmax,ifirst
	data ifirst/.true./
	data epsilon/0.001/
c
c	ICODE = 0       no geometry
c             = 1	rectangular geometry
c	      = 2	cylindrical geometry
c       IFLAG = 0       the track segment is all inside inner cylinder
c             = 1       the track segment crosses the inner/outer wall
c             = 2       the track segment crosses the outer wall boundary
c             = 3       the track segment is all outside outer cylinder
c
	if(ifirst)then
	   if(icode.eq.1)then
	      write(6,*)'Geometry of target inconsisitent with code'
	      call exit ()
	   elseif(icode.eq.2)then
	      radius1 =  size1  * 0.5
	      radius2 =  size2  * 0.5
	      zmin    =  size3
	      zmax    =  size4
	      ifirst  = .false.
c	      write(6,*)radius1,radius2,zmin,zmax
	   elseif(icode.eq.0)then
	      return
	   else
	      write(6,*)'Target geometry: huh?'
	      call exit()
	   endif
	endif
c       
c	Is the starting point for this track segment inside the cylinder?
c
	istart = -1                         !flags where segment starts
	iend   = -1
	iflag  = -1
	rr0 = sqrt(r0(1)**2. + r0(2)**2.)
	if(rr0.lt.radius2 .and.
     1               r0(3).ge.zmin 
     1         .and. r0(3).le.zmax )then
	   if(rr0.lt.(radius1-epsilon) .and.
     1  	      r0(3).ge.(zmin+epsilon) 
     1          .and. r0(3).le.(zmax-epsilon) )then
	      istart = 0
	   else
	      istart = 1
	   endif
	else
	   istart = 2
	endif
c	write(6,*)'>   istart,rr0,radius1,radius2,zmin,zmax'
c	write(6,*)'>  ',istart,rr0,radius1,radius2,zmin,zmax
c
c	Check if the ending point is outside cylinder
c
	rr1 = sqrt(r1(1)**2. + r1(2)**2.)
	if(rr1.lt.radius2 .and.
     1                r1(3).ge.zmin
     1        .and.  r1(3).le.zmax )then
	   if(rr1.lt.(radius1-epsilon) .and.
     1               r1(3).ge.(zmin+epsilon)
     1         .and. r1(3).le.(zmax-epsilon) )then
	      iend  = 0
	   else
	      iend  = 1
	   endif
	else	
	   iend  = 2
	endif
	if(istart.eq.-1.or.iend.eq.-1)then
	   write(6,*)'Deep do-do',r0,r1,rr0,rr1
	end if
c	write(6,*)'>> iend,rr1'
c	write(6,*)'>> ',iend,rr1
	if(istart.eq.0 .and. iend.eq.0)then
	   iflag = 0
	   return
	elseif (istart.eq.0 .and. iend.eq.1)then
	   iflag = 1
	   call linecircle(r0(1),r0(2),r1(1),r1(2),radius1,xx,yy)
	elseif (istart.eq.0 .and. iend.eq.2)then
	   iflag = 1
	   call linecircle(r0(1),r0(2),r1(1),r1(2),radius1,xx,yy)
	elseif (istart.eq.1 .and. iend.eq.2)then
	   iflag = 2
	   call linecircle(r0(1),r0(2),r1(1),r1(2),radius2,xx,yy)
	elseif (istart.eq.2 .and. iend.eq.2)then
	   iflag = 3          !the whole segment is all outside
	   return
	elseif (istart.eq.1 .and. iend.eq.1)then
	   iflag = 0
	   return
	elseif (istart.eq.2 .and. iend.eq.0)then
	   iflag = 3  !rare backward-traveling particle; ignore
	   return
	elseif (istart.eq.1 .and. iend.eq.0)then
	   iflag = 3  !rare backward-traveling particle; ignore
	   return
	elseif (istart.eq.2 .and. iend.eq.1)then
	   iflag = 3  !rare backward-traveling particle; ignore
	   return
	endif
	if(iflag.eq.-1)then
	   write(6,*)'Help!,istart,iend,r0,r1'
	   write(6,*)'Help!',istart,iend,r0,r1
	endif
	zz = r0(3) + direc(3)/direc(1) * (xx-r0(1))
c	write(6,10)'>>> r0,r1,istart,iend,iflag,xx,yy,zz'
c	write(6,10)'>>>',r0,r1,istart,iend,iflag,xx,yy,zz
c 10	format(1x,a,3f6.2,3f6.2,3i3,3f6.2)
	r1(1) = xx 
	r1(2) = yy 
	r1(3) = zz
	if(r1(3).gt.zmax)then
	   ratio = (zmax - r0(3)) / (r1(3) - r0(3))
	   r1(3) = zmax 
	   r1(1) = r0(1) + ratio * (r1(1) - r0(1))
	   r1(2) = r0(2) + ratio * (r1(2) - r0(2))
	else if (r1(3).lt.zmin)then
	   ratio = (zmin - r0(3)) / (r1(3) - r0(3))
	   r1(3) = zmin 
	   r1(1) = r0(1) + ratio * (r1(1) - r0(1))
	   r1(2) = r0(2) + ratio * (r1(2) - r0(2))
	end if
	deltad = sqrt((r1(1)-r0(1))**2.+(r1(2)-r0(2))**2.+(r1(3)-r0(3))**2.)
	return
	end
	
c
c     Find intersection of a line with a circle
c     Do it using a parametric equation for the line in "t"
c
	subroutine linecircle(x0,y0,x1,y1,radius,xx,yy)
c       
	theta = atan2(y1-y0,x1-x0)
	cost = cos(theta)
	sint = sin(theta)
c       a    = 1.
c       twoa = 2.
	b    = 2.*(y0*sint+x0*cost)
	c    = x0*x0 + y0*y0 - radius*radius
	rad  = sqrt(b*b - 4.*c)
	tplus  = 0.5 * (-b + rad)
c       tminus = (-b - rad) / twoa
c       write(6,*)tplus,tminus
	xx = x0 + tplus*cost
	yy = y0 + tplus*sint
	return
	end

	

	


