
	subroutine init_parms(filename)
C
	implicit none
C
	include 'constants.inc'
C
	integer*4   lenocc
	external    lenocc
C
	integer*4   db_maxsize,db_nfields
	parameter   (db_maxsize=80)
	parameter   (db_nfields=1)

	integer*4   j, ident, ind, avail, iok
	real*8	    db_field(db_maxsize,db_nfields)
	logical     success
	character*20   filename
	character*80   string
	character*20   ftitle
	character*10   phbstr
	integer*4      ipos,jpos
C
	save
C
	filename = filename(1:lenocc(filename))
	open(unit=1, file=filename, status='old',iostat=iok)
	if(iok.ne.0) then
	   write(6,*) ' That database doesn''t exist here! Try again ...'
	   print *,filename
	   stop 'Exiting ...'
	endif	  

C Read in the dbase fields

C ... watch for the ident marker
10	read(1,'(a)',end=70,err=60) string
	read(string(41:60),*) ftitle
	ftitle = ftitle(1:lenocc(ftitle))
	read(string(1:6),1010) ident
 1010	format(I7)
	ind = 32768+ident
C
	if(ind.gt.0.and.ind.le.db_nfields) goto 30
	write(6,'(1x,a)') 'WARNING: One of the dbase fields is the wrong size'
	goto 10

C ... read in field and go back for the next one
 30	read(string(8:9),1011) avail
 1011	format(I2)
C 	write(6,'(1x,a,a,a)') '\n Reading in DBASE field ... ',ftitle,avail

	do j=1,avail
	   read(1,*,end=50,err=60) db_field(j,ind)
C	   write(6,*) 'Field',j,ind,db_field(j,ind)
	enddo
C
	goto 10

C ... oops
50	write(6,'(1x,a,i3)')
	1	'ERROR: read_params found EOF in the middle of field #', ind
	close(1)
	return
60	write(6,'(1x,a)') 'ERROR: Read error in read_params' 
	close(1)
	return
70	if (ind.ne.db_nfields) then
	  write(6,'(1x,a,i2,a,i2)') 'ERROR: read_params found ',ind,
	1	' fields instead of ',db_nfields
	  close(1)
	  return

C ... success!
	else
	  success=.true.
	  close(1)
	endif

C Grab what we want from the information

C ... dbase field 
	ind=1
	omega2max = db_field(1,ind)
	ncall = int(db_field(2,ind))
	itmx(1) = int(db_field(3,ind))
	acc(1) = db_field(4,ind)
	itmx(2) = int(db_field(5,ind))
	acc(2) = db_field(6,ind)
	itmx(3) = int(db_field(7,ind))
	acc(3) = db_field(8,ind)
	itmx(4) = int(db_field(9,ind))
	acc(4) = db_field(10,ind)
	mxevnt = int(db_field(11,ind))
	if (db_field(12,ind).EQ.1.0) then
	   nton = .TRUE.
	else
	   nton = .FALSE.
	endif
C
	return
C
	END


