	subroutine NtupleFill(iun)

	implicit none
C
        include 'hbook.inc'
	include 'kinem.cmn'
	include 'constants.inc'
C
	integer	iun
	real*4	ntu(nentries)
C
	integer*4 j 
C
	save
C
	j = 0
C
	j = j+1
	ntu(j) = iEbin
C	j = j+1
C	ntu(j) = ievent
	j = j+1
	ntu(j) = sngl(Egamma)
	j = j+1
	ntu(j) = sngl(Eprime)
	j = j+1
	ntu(j) = sngl(costheta)
	j = j+1
	ntu(j) = sngl(phi)
	j = j+1
	ntu(j) = sngl(Egrad)
	j = j+1
	ntu(j) = sngl(cthgrad)
	j = j+1
	ntu(j) = sngl(phir)
	j = j+1
	ntu(j) = sngl(phirel)
	j = j+1
	ntu(j) = sngl(cthrel)
	j = j+1
	ntu(j) = sngl(Eelec)
	j = j+1
	ntu(j) = sngl(phie)
	j = j+1
	ntu(j) = sngl(costhe)
	j = j+1
	ntu(j) = sngl(xgamma)
	j = j+1
	ntu(j) = sngl(ygamma)
	j = j+1
	ntu(j) = sngl(xelec)
	j = j+1
	ntu(j) = sngl(yelec)
	j = j+1
	ntu(j) = sngl(xgrad)
	j = j+1
	ntu(j) = sngl(ygrad)
C	j = j+1
C	ntu(j) = sngl(dx)
C	j = j+1
C	ntu(j) = sngl(dy)
C	j = j+1
C	ntu(j) = sngl(rmin)
	j = j+1
	ntu(j) = sngl(corr)
	j = j+1
	ntu(j) = sngl(hardcorr)
	j = j+1
	ntu(j) = sngl(vsoftcorr)
C	j = j+1
C	ntu(j) = sngl(Pxtot)
C	j = j+1
C	ntu(j) = sngl(Pytot)
C	j = j+1
C	ntu(j) = sngl(Pztot)
C	j = j+1
C	ntu(j) = sngl(Pg1(1))
C	j = j+1
C	ntu(j) = sngl(Pg1(2))
C	j = j+1
C	ntu(j) = sngl(Pg1(3))
C	j = j+1
C	ntu(j) = sngl(Pg1(4))
C	j = j+1
C	ntu(j) = sngl(Pg2(1))
C	j = j+1
C	ntu(j) = sngl(Pg2(2))
C	j = j+1
C	ntu(j) = sngl(Pg2(3))
C	j = j+1
C	ntu(j) = sngl(Pg2(4))
C	j = j+1
C	ntu(j) = sngl(Pe1(1))
C	j = j+1
C	ntu(j) = sngl(Pe1(2))
C	j = j+1
C	ntu(j) = sngl(Pe1(3))
C	j = j+1
C	ntu(j) = sngl(Pe1(4))
	j = j+1
	ntu(j) = cutg
	j = j+1
	ntu(j) = cute
	j = j+1
	ntu(j) = cutr
C	j = j+1
C	ntu(j) = e1cut
C	j = j+1
C	ntu(j) = e2cut
C	j = j+1
C	ntu(j) = p1cut
!
        call HFN(NtupleID,ntu)        
	return
	end
	
