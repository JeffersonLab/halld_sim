c******************************************************************************
c	"CARROT.FOR"
c	CARROT does a rotation of the Cartesian vector VEC by ROTANG
c	radians about the arbitrary Cartesian vector ROT.  The 
c	rotation is clockwise when viewed along the rotation axis.
c
c	R.A.Schumacher, S.I.N., Dec. 86
c
c       9-97 Unix version - got rid of real*4
c       - tweaked cross product routine for speed
c	6-28-17 - make vectors of length 5 for compatibility with mc_gen changes
c                 even though we are doing simple 3-vector operations
c	
	subroutine carrot(vec,rot,rotang)
	real vec(5),rot(5),rhat(5),vecpara(5),vecperp(5)
	real vcross(5),vecperprot(5)
c
	rr = sqrt(rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3))
	if(rr.eq.0)then
		rotang = 0.0
		return
	end if
	rhat(1) = rot(1)/rr		!direction cosines of rotation axis
	rhat(2) = rot(2)/rr
	rhat(3) = rot(3)/rr
c
c	Decompose VEC into components parallel and perpendicular to ROT
c
	vecdot     = vec(1)*rhat(1) + vec(2)*rhat(2) + vec(3)*rhat(3)
	vecpara(1) = vecdot * rhat(1)	!projection of VEC along rotation axis
	vecpara(2) = vecdot * rhat(2)
	vecpara(3) = vecdot * rhat(3)
	vecperp(1) = vec(1) - vecpara(1)!projection of VEC perpendicular to 
	vecperp(2) = vec(2) - vecpara(2)!the rotation axis
	vecperp(3) = vec(3) - vecpara(3)
c
c	The expression for VPERPROT is derived from the cross product
c	of VPERP and VPERPROT, and the "BAC-CAB" identity.
c
	call crossprod(vecperp,rhat,vcross,angle)
	sinrot = sin(rotang)		!rotation angle about ROT
	cosrot = cos(rotang)
	vecperprot(1) = vecperp(1)*cosrot - vcross(1)*sinrot
	vecperprot(2) = vecperp(2)*cosrot - vcross(2)*sinrot
	vecperprot(3) = vecperp(3)*cosrot - vcross(3)*sinrot
c
	vec(1) = vecperprot(1) + vecpara(1)	!Resultant VEC
	vec(2) = vecperprot(2) + vecpara(2)
	vec(3) = vecperprot(3) + vecpara(3)
	continue
	return
	end

c
c	Cross product subroutine.
c	A x B = C; ANGLE is the angle between A and B in radians.
c
	subroutine crossprod(a,b,c,angle)
	real a(5),b(5),c(5)
c
	c(1) =   b(3)*a(2) - b(2)*a(3)
	c(2) = -(b(3)*a(1) - b(1)*a(3))
	c(3) =   b(2)*a(1) - b(1)*a(2)
	c(4) = 0.
	c(5) = 0.
	aa = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
	bb = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
	cc = c(1)*c(1) + c(2)*c(2) + c(3)*c(3)
	dd = aa*bb
	if(dd.ne.0.0)then
		sinang = sqrt(cc/dd)
		if(sinang.gt.1.0)sinang=1.0
	else
		sinang = 0.0
	end if
c	write(6,*)'Cross prod:',a
c	write(6,*)'Cross prod:',b
c	write(6,*)'Cross prod:',c
c	write(6,*)'Cross prod:',aa,bb,cc,sinang
	angle  = asin(sinang)
	return
	end

c
c	Test program
c
c	real oldvec(3),vec(3),rot(3)
c10	read(5,*)vec,rot,ang
c	rotang = ang*3.1415926535/180.
c	oldvec(1) = vec(1)
c	oldvec(2) = vec(2)
c	oldvec(3) = vec(3)
c	call carrot(vec,rot,rotang)
c	write(6,*)'OLDVEC',oldvec
c	write(6,*)'   VEC',vec
c	write(6,*)'   ROT',rot
c	write(6,*)' ANGLE',ang
c	goto 10
c	end
