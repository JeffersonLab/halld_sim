c*******************************************************************************
c       Relativistic lineshape with phase space and Breit-Wigner: random numbers.
c       Returns a mass value drawn from a relativistic Breit-Wigner lineshape
c       with a Gamma factor that depends on the mass.
c
c	This routine includes the initial state phase space factor as well, which
c	distribguishes it from its parent routine "relativisticbwan.f".
c
c	It is a complete rewrite-however, since we can't use the integration methods
c	when the probablility density function is different event by event.
c	
c       Inputs:
c       rmnot    - central mass
c       gammanot - nominal width of the state
c       rm1      - first  decay product mass (needed to compute momentum)
c       rm2      - second decay product mass
c	rm3      - the particle recoiling against the resonance with mass rmnot
c       lvalue   - integer value specifying orbital L value of decay (0,1,2...)
c	ppbeam   - beam energy in lab
c
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
c       01-21-11 - modified to give a momentum-independent lineshape in cases when
c                  the requested decay products are more massive than the centroid
c                  of the parent particle
c	08-29-17 - Introduce the phase space factor for the recoiling particle in the CM
c                 frame, e.g. the recoiling target proton in the case of pbar-p creation
c	09-02-17 - this is a rather complete re-write to adopt the "dart-throwing" monte-carlo method
c	          rather than the integration method... on account of every call being for
c	          a different value of W, so we cannot pre-compute the entire line shape.
c       09-22-17 - If the width is "-9999", return just the value drawn from phase-space
c       03-06-19 - Trap possible failure in computing RMMAX
c       07-17-19 - Improve comments written to screen
c       11-08-19 - turn off some rare and not helpful error messages
c	10-12-20 - correct /monte/ common block and declare "external ran"
c       
	subroutine relativisticlineshape(rmnot,gammanot,rm1,rm2,lvalue,
     1                             	rm3,ww,bwmass,ifail)
c       Inputs:
c	RMNOT     - Mass of centroid of relativisitc Breit Wigner being generated
c	GAMMANOT  - Nominal width of the resonance, before q-dependence inclusion
c	RM1       - Mass of the first decay product
c	RM2       - Mass of the second decay product
c	LVALUE    - Orbital angular momentum of the decaying pair - for Gamma(q)
c	RM3       - Mass of recoiling particle in 3-body final state - for phase space calc.
c	WW        - Total available energy
c
c	Outputs:
c	BWMASS    - Monte Carlo selection of mass 
c	IFAIL     - error flag
c
	external ran
	logical ifirst
	data ifirst/.true./
	common /monte/ iseed,xflat,xnormal,xbw,ixcount
	common /lineshapelocal/ifirst,bwmaxval,rmnotsave,scale
c
	ifail = 0
	if(ifirst.eqv..true.)then
	   ifirst    = .false.
c	   write(6,*)' '
	   write(6,*)'First call to Relativisitic Breit-Wigner with phase space:'
	   write(6,10)rmnot,gammanot,rm1,rm2,rm3,ww,lvalue
 10	   format(1x,
     1     'Central resonance mass          ',f8.3,' GeV'/,1x,
     1     'Nominal width                   ',f8.3,/,1x,
     1     'Mass of first decay product     ',f8.3,/,1x,
     1     'Mass of second decay product    ',f8.3,/,1x,
     1     'Mass of recoiling particle      ',f8.3,/,1x,
     1     'Total available energy (W)      ',f8.3,/,1x,
     1     'Orbital angluar momentum (L)    ',i8)
c	   write(6,*)' '
	   rmnotsave = rmnot
	   bwmaxval = 1000.
	   scale = 1000.
c	   call exit()
	endif
c
c       Establish the limits of the mass generation to be between
c       4 gamma below to 20 gamma above the nominal centroid,
c       as long as these masses are legal.
c
	rmmin     = max(rmnot -  4.*gammanot,rm1+rm2)
	rmmax     = min(ww-rm3, rmnot + 20.*gammanot)
	drm       = rmmax-rmmin
	if(drm.lt.0)then
	   ifail = 1
	   bwmass = rmnot
	   return
	endif
c       
c       Pick a random mass value to test against the BW line shape
c       The limits of allowed masses are already established
c
	icount = 0
 100	icount = icount + 1
	rm     = rmmin + ran(iseed)*drm
c
c       Compute the momentum between the created resonance, of mass
c       RM, and the recoiling particle RM3
c
	rad = (ww**2.-(rm**2.+rm3**2.))**2. - (2.*rm*rm3)**2.
	if(rad.ge.0.0)then
	   pcm = sqrt(rad)/(2.*ww)
	else
c	   write(6,*)'P bad...',ww,rmnot,rm3,rm1,rm2,rm
c	   write(6,*)rmmin,rmmax,rad
	   pcm = 0.0
	   ifail = 3
	   bwmass = rmnot
	   return
c           call exit()
	endif
c
c       Scale the first part of the phase space calculation using the
c       CM momentum evaluated at the threshold of the reaction
c       (It makes no difference which mass is used when computing the ratio.)
c
	rad = (ww**2.-(rmmin**2.+rm3**2.))**2. - (2.*rmmin*rm3)**2.
c	rad = (ww**2.-(rmnot**2.+rm3**2.))**2. - (2.*rmnot*rm3)**2.
	if(rad.ge.0.0)then
	   pcmnot = sqrt(rad)/(2.*ww)
	else
	   write(6,*)'P0 bad...',ww,rmnot,rm3,rm1,rm2,rm
	   write(6,*)rmmin,rmmax,rad
	   ifail = 4
	   bwmass = rmnot
	   return
c           call exit()
	endif
	pratio = pcm/pcmnot
c
c       Compute the momentum of the final state into which the resonance decays
c
	rad = (rm**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	if(rad.ge.0.0)then
	   qcm     = sqrt(rad)/(2.*rm)
	else
	   write(6,*)'Q bad...',ww,rmnot,rm3,rm1,rm2,rm
	   ifail = 2
	   bwmass = rmnot
	   return
c	   call exit()
	endif
	rad = (rmnot**2.-(rm1**2.+rm2**2.))**2. - (2.*rm1*rm2)**2.
	if(rad.ge.0.0)then
	   qcmnot = sqrt(rad)/(2.*rmnot)
	else
	   write(6,*)'Q0 bad...',ww,rmnot,rm3,rm1,rm2,rm
	   write(6,*)rmmin,rmmax,rad
	   ifail = 2
	   bwmass = rmnot
	   return
c	   call exit()
	endif
	qratio  = qcm/qcmnot
c       
c       Compute the value of the line shape at this selected mass RM
c
	rmratio = rmnot/rm
	gamma   = gammanot * rmratio*qratio**(2.*lvalue+1.)
	phasespace = qratio * pratio
	bwvalue = scale * phasespace
	if(gammanot.lt.9998)then  !special flag to return just phase space
	   bwvalue = bwvalue*(rm*rmnot*gammanot*qratio**(2.*lvalue))/
     1	          ((rm**2. - rmnot**2.)**2. + (rmnot*gamma)**2.)
	endif
c       
c       Throw the dart...
c
	ytest = ran(iseed) * bwmaxval
c	write(6,*)'Generating a line shape: ',icount,rm,ytest,bwvalue,bwmaxval
	if(bwvalue.gt.bwmaxval)then
	   bwmaxval = bwvalue
c	   write(6,*)'New maximal value in BW lineshape:',bwmaxval,qratio,pratio
	endif
	if(ytest.gt.bwvalue)then
c	   write(6,*)'Go back'
	   if(icount.lt.1000)then
	      goto 100
	   else
	      write(6,*)'Inefficient BW shape',icount,rmnot,rm1,rm2,rm3,qcm,pcm
c	      call exit()
	   endif
	endif
	bwmass = rm
c	write(6,*)'Relativistic BW lineshape:',bwmass
	return
	end

