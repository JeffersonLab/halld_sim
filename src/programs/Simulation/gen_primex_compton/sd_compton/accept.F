	subroutine accept
C
	implicit none
C
	include 'kinem.cmn'
	include 'constants.inc'
C
	save
C
	e1cut = 0
	e2cut = 0
	p1cut = 0
C
	cutg = 0
	cutr = 0
	cute = 0
C
	if (Eprime.GT.HyCal_thresh) then
	   e1cut = 1
	   if ((abs(xgamma).GT.genlim(1).OR.abs(ygamma).GT.genlim(3)).AND.
	1	(abs(xgamma).LT.genlim(2).AND.abs(ygamma).LT.genlim(4))) cutg = 1
	endif
C
	if (Egrad.GT.HyCal_thresh) then
	   e2cut = 1
	   if ((abs(xgrad).GT.genlim(1).OR.abs(ygrad).GT.genlim(3)).AND.
	1	(abs(xgrad).LT.genlim(2).AND.abs(ygrad).LT.genlim(4))) cutr = 1
	endif
C
	if (Eelec.GT.HyCal_thresh) then
	   p1cut = 1
	   if ((abs(xelec).GT.genlim(1).OR.abs(yelec).GT.genlim(3)).AND.
	1	(abs(xelec).LT.genlim(2).AND.abs(yelec).LT.genlim(4))) cute = 1
	endif
C
	RETURN
	END
