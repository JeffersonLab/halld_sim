c*******************************************************************************
c       Relativistic Breit-Wigner random numbers.
c       Returns a mass value drawn from a relativistic Breit-Wigner lineshape
c       with a Gamma factor that depends on the mass.  The first call 
c       sets up a mapping table to choose the mass within -5 +30 Gamma(m_not)
c       values of the centroid mass.  Later calls use the established mapping.
c       That is, this routine uses the "transformation" method decribed, 
c       e.g.,  in Bevington page 82).  
c	
c       Inputs:
c       rmnot    - central mass
c       gammanot - nominal width of the state
c       rm1      - first  decay product mass (needed to compute momentum)
c       rm2      - second decay product mass
c       lvalue   - integer value specifying orbital L value of decay (0,1,2...)
c       Outputs:
c       bwmass   - randomized mass value generated from the function 
c                  being modeled
c       
c	R.A.Schumacher, CMU, 11-2008.
c
c       11-18-08 - started with BWRAN code and made major changes
c       11-24-08 - widen the range over which the lineshape is computed.  Make
c                  the range asymmetric, reflecting the fact that the long tail
c                  is on the high side of the distribution
c        1-21-11 - modified to give a momentum-idependent lineshape in cases when
c                  the requested decay products are more massive than the centroid
c                  of the parent particle
c	 8-10-17 - turned this into TWO identical subroutines so that I can call two
c                  separate instances within a single reaction description
c         	   (Lambda anti-Lambda case)  Only the common blocks differ
c       10-12-20 - update /monte/ common block and declare "external ran"
	subroutine relativisticbwran1(rmnot,gammanot,rm1,rm2,lvalue,bwmass)
	external ran
	parameter (nbin=1000)
	real bwint(0:nbin)
	logical ifirst,ignoremasses
	data ifirst/.true./
	common /monte/ iseed,xflat,xnormal,xbw,ixcount
	common /rbwlocal1/ifirst,ignoremasses,
     1                   bwint,rmmin,drm,rmnotsave,rm1save,rm2save
	if(ifirst.eqv..true.)then
	   ifirst    = .false.
c	   write(6,*)' '
	   write(6,*)'Initializing Relativisitic Breit-Wigner', rm1, rm2
	   write(6,*)'Orbital angular momentum is ',lvalue
	   write(6,*)' '
	   rmnotsave = rmnot
	   rm1save   = rm1
	   rm2save   = rm2
	   rmmin     = rmnot -  4.*gammanot
	   rmmax     = rmnot + 20.*gammanot
	   drm       = (rmmax-rmmin)/nbin
	   bwint(0)  = 0.0
c
c          'Integrate' the function and store the partial sums
c
	   gamma = gammanot  !plain-old mass-independent Breit-Wigner gamma value
	   rad = (rmnot**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	   if(rad.ge.0.0)then
	      pcmnot = sqrt(rad)/(2.*rmnot)
	      ignoremasses = .false.
	   else
	      write(6,*)'Central mass value is below threshold',rmnot,rm1,rm2
	      pcmnot = 1.0
	      write(6,*)'******************************************************'
	      write(6,*)'The calculation will use NON-relativistic Breit-Wigner'
	      write(6,*)'        ...we hope you know what you are doing.'
	      write(6,*)'******************************************************'
	      ignoremasses = .true.
c	      call exit()
	   endif
c
	   do i= 1,nbin
	      rm = rmmin + (i-1) * drm
	      rad = (rm**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	      if(rad.ge.0.0 .and. .not.ignoremasses .and. rm.gt.0.0)then
		 pcm     = sqrt(rad)/(2.*rmnot)
		 gamma   = gammanot * (pcm/pcmnot)**(2.*lvalue+1.)
		 bwvalue = 2.*rm*rmnot*gamma/
     1	               ((rm**2. - rmnot**2.)**2. + (rmnot*gamma)**2.)
	      else
c		 write(6,*)'Below threshold',i,rmnot,rm,rm1,rm2
		 gamma   = gammanot
		 pcm = 0.0
		 bwvalue = 0.0
	      endif
 	      bwint(i) = bwint(i-1) + bwvalue
c	      write(6,*)'Values:        ',i,rmnot,rm,rm1,rm2,pcm,gamma,bwvalue
 	   enddo
c
c          Normalize the integral
c
	   fact=bwint(nbin)  !last bin contains the unnormalized total
	   do i = 1,nbin
	      bwint(i) = bwint(i)/fact
	   enddo
	end if 
c	call exit()
c       
c       If this routine was called with "new" masses, we are not prepared to 
c       deal with it.
c
	if(rmnot.ne.rmnotsave .or. rm1.ne.rm1save .or. rm2.ne.rm2save)then
	   write(6,*)'Mismatch of masses encountered in RelatvisiticBWran'
	   write(6,*)rmnot,rm1,rm2,rmnotsave,rm1save,rm2save
	   call exit()
	endif
c
	icount = 0
 10	icount = icount+1
	test = ran(iseed)
	do i=1,nbin
c	   write(6,*)i,bwint(i),test,drm,rmmin
	   if(bwint(i).ge.test)then
	      bwmass = rmmin + (i-1) * drm
c
c	      We blur the result a bit, within one bin
c             to avoid some binning problems
c
	      bwmass = bwmass + drm*(ran(iseed)-0.5)
	      goto 40
	   end if
	enddo
	write(6,*)'Bloody Murder',i,bwmass,test
	call exit()
 40	if(bwmass.lt.(rm1+rm2) .and. .not.ignoremasses)then
	   if(icount.lt.20)then
	      goto 10
	   else
	      write(6,*)'Below threshold mis-match in RelativisticBWran'
	      write(6,*)bwmass,rm1,rm2
	      call exit()
	   endif
	endif
c	write(6,*)'Relativistic Breit-Wigner mass:',bwmass,rmnot,rm1,rm2,icount
	return 
	end

C       Second identical instance of this subroutine (excpet for separate common block)
C
	subroutine relativisticbwran2(rmnot,gammanot,rm1,rm2,lvalue,bwmass)
	external ran
	parameter (nbin=1000)
	real bwint(0:nbin)
	logical ifirst,ignoremasses
	data ifirst/.true./
	common /monte/ iseed,xflat,xnormal,xbw,ixcount
	common /rbwlocal2/ifirst,ignoremasses,
     1                   bwint,rmmin,drm,rmnotsave,rm1save,rm2save
	if(ifirst.eqv..true.)then
	   ifirst    = .false.
c	   write(6,*)' '
	   write(6,*)'Initializing Relativisitic Breit-Wigner', rm1, rm2
	   write(6,*)'Orbital angular momentum is ',lvalue
	   write(6,*)' '
	   rmnotsave = rmnot
	   rm1save   = rm1
	   rm2save   = rm2
	   rmmin     = rmnot -  4.*gammanot
	   rmmax     = rmnot + 20.*gammanot
	   drm       = (rmmax-rmmin)/nbin
	   bwint(0)  = 0.0
c
c          'Integrate' the function and store the partial sums
c
	   gamma = gammanot  !plain-old mass-independent Breit-Wigner gamma value
	   rad = (rmnot**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	   if(rad.ge.0.0)then
	      pcmnot = sqrt(rad)/(2.*rmnot)
	      ignoremasses = .false.
	   else
	      write(6,*)'Central mass value is below threshold',rmnot,rm1,rm2
	      pcmnot = 1.0
	      write(6,*)'******************************************************'
	      write(6,*)'The calculation will use NON-relativistic Breit-Wigner'
	      write(6,*)'        ...we hope you know what you are doing.'
	      write(6,*)'******************************************************'
	      ignoremasses = .true.
c	      call exit()
	   endif
c
	   do i= 1,nbin
	      rm = rmmin + (i-1) * drm
	      rad = (rm**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	      if(rad.ge.0.0 .and. .not.ignoremasses .and. rm.gt.0.0)then
		 pcm     = sqrt(rad)/(2.*rmnot)
		 gamma   = gammanot * (pcm/pcmnot)**(2.*lvalue+1.)
		 bwvalue = 2.*rm*rmnot*gamma/
     1	               ((rm**2. - rmnot**2.)**2. + (rmnot*gamma)**2.)
	      else
c		 write(6,*)'Below threshold',i,rmnot,rm,rm1,rm2
		 gamma   = gammanot
		 pcm = 0.0
		 bwvalue = 0.0
	      endif
 	      bwint(i) = bwint(i-1) + bwvalue
c	      write(6,*)'Values:        ',i,rmnot,rm,rm1,rm2,pcm,gamma,bwvalue
 	   enddo
c
c          Normalize the integral
c
	   fact=bwint(nbin)  !last bin contains the unnormalized total
	   do i = 1,nbin
	      bwint(i) = bwint(i)/fact
	   enddo
	end if 
c	call exit()
c       
c       If this routine was called with "new" masses, we are not prepared to 
c       deal with it.
c
	if(rmnot.ne.rmnotsave .or. rm1.ne.rm1save .or. rm2.ne.rm2save)then
	   write(6,*)'Mismatch of masses encountered in RelatvisiticBWran'
	   write(6,*)rmnot,rm1,rm2,rmnotsave,rm1save,rm2save
	   call exit()
	endif
c
	icount = 0
 10	icount = icount+1
	test = ran(iseed)
	do i=1,nbin
c	   write(6,*)i,bwint(i),test,drm,rmmin
	   if(bwint(i).ge.test)then
	      bwmass = rmmin + (i-1) * drm
c
c	      We blur the result a bit, within one bin
c             to avoid some binning problems
c
	      bwmass = bwmass + drm*(ran(iseed)-0.5)
	      goto 40
	   end if
	enddo
	write(6,*)'Bloody Murder',i,bwmass,test
	call exit()
 40	if(bwmass.lt.(rm1+rm2) .and. .not.ignoremasses)then
	   if(icount.lt.20)then
	      goto 10
	   else
	      write(6,*)'Below threshold mis-match in RelativisticBWran'
	      write(6,*)bwmass,rm1,rm2
	      call exit()
	   endif
	endif
c	write(6,*)'Relativistic Breit-Wigner mass:',bwmass,rmnot,rm1,rm2,icount
	return 
	end


