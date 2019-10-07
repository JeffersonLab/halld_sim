c*******************************************************************************
c	Breit-Wigner distributed random numbers.
c	BW is distributed about zero in the range -5 to +5 according to
c	the probability density function 1/(BW**2 + 1), and normalized.
c	Take the BW returned by this subroutine and turn it into an energy
c	(or mass) by computing Mass = (Gamma/2)*BW + Epole, where 
c	Gamma/2 is the width of the resonance, and Epole is the central
c	energy of the resonance.
c
c	R.A.Schumacher, CMU, 12-92
c
c       8-18-97  - UNIX version
c       3-18-99  - LINUX version
c       11-21-05 - improve code, broadened range from +-5 to +-10.
c       1-17-11  - shift  range to be -5  to +15
c       5-6-11   - expand range to be -15 to +15
c       5-11-11  - expand range to be +30 to -30
c
	subroutine bwran(bw)
	parameter (nbin=1000)
	real bwi(0:nbin)
	logical ifirst
	data ifirst/.true./
	common /monte/ iseed,xflat,xnormal,xbw
	common /bwlocal/bwi,xmin,dx !needed for Linux fortran
	if(ifirst.eqv..true.)then
	   ifirst = .false.
	   dx     = 60.0/nbin
	   xmin   = -30.0
	   bwi(0) = 0.0
	   do 10 i= 1,nbin
	      x = xmin + (i-1) * dx
 	      bwi(i) = bwi(i-1) + 1.0/(x*x + 1.0)
 10	   continue
	   fact=bwi(nbin)
	   do 20 i= 1,nbin
	      bwi(i) = bwi(i)/fact
 20	   continue
	end if 
	test = ran(iseed)
	do 30 i=1,nbin
	   if(bwi(i).ge.test)then
	      bw = xmin + (i-1) * dx
c
c	      We blur the result a bit, within one bin
c             to avoid some binning problems
c
	      bw = bw + dx*(ran(iseed)-1.0)
	      return
	   end if
 30	continue
	xbw = bw		!for diagnostic plotting
	write(6,*)'Bloody Murder',i,bw,test
	call exit()
	end
