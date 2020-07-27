      subroutine NtupleInit(filename)      
      implicit none
      save
*
      include  'hbook.inc'
      include 'constants.inc'
*
      character*80 filename,directory
      character*16 NtupleTag(nentries),name,title
      integer*4 m,io,recl,bank,id,NtupleSize,status
C      parameter(recl = 1024)
      parameter(recl = 8190)
      parameter(bank = 8000)
      parameter(io = 27)
      parameter(name = 'RCNtuple')
      parameter(title = 'COMP NTUPLE')
C
      integer*4    lenocc
      external     lenocc
*
      NtupleID = defaultID
      id = NtupleID
      NtupleIO = io
      NtupleName = name
*
      call HCDIR(directory,'R') !CERNLIB read current directory
      call HROPEN(io,name,filename,'N',recl,status)  !CERNLIB
*                                       !directory set to "//TUPLE"
      if (status.ne.0)then
         write(6,*)'HROPEN error: istat=',status
         stop
      else
         write(6,*)' '
         write(6,*)'Opened file: ',filename(1:lenocc(filename))
      endif
      m = 0
*      
      m = m+1
      NtupleTag(m) = 'Ebin'
C      m = m+1
C      NtupleTag(m) = 'event'
      m = m+1
      NtupleTag(m) = 'Ebeam'
      m = m+1
      NtupleTag(m) = 'Egamma2'
      m = m+1
      NtupleTag(m) = 'ctheta2'
      m = m+1
      NtupleTag(m) = 'phi2'
      m = m+1
      NtupleTag(m) = 'Egamma1'
      m = m+1
      NtupleTag(m) = 'ctheta1'
      m = m+1
      NtupleTag(m) = 'phi1'
      m = m+1
      NtupleTag(m) = 'phirel'
      m = m+1
      NtupleTag(m) = 'cthrel'
      m = m+1
      NtupleTag(m) = 'Ee'
      m = m+1
      NtupleTag(m) = 'phie'
      m = m+1
      NtupleTag(m) = 'costhe'
      m = m+1
      NtupleTag(m) = 'xgamma2'
      m = m+1
      NtupleTag(m) = 'ygamma2'
      m = m+1
      NtupleTag(m) = 'xelec'
      m = m+1
      NtupleTag(m) = 'yelec'
      m = m+1
      NtupleTag(m) = 'xgamma1'
      m = m+1
      NtupleTag(m) = 'ygamma1'
C      m = m+1
C      NtupleTag(m) = 'dx'
C      m = m+1
C      NtupleTag(m) = 'dy'
C      m = m+1
C      NtupleTag(m) = 'rmin'
      m = m+1
      NtupleTag(m) = 'cswitch'
      m = m+1
      NtupleTag(m) = 'hardcorr'
      m = m+1
      NtupleTag(m) = 'vscsec'
C      m = m+1
C      NtupleTag(m) = 'Pxtot'
C      m = m+1
C      NtupleTag(m) = 'Pytot'
C      m = m+1
C      NtupleTag(m) = 'Pztot'
C      m = m+1
C      NtupleTag(m) = 'Pg11'
C      m = m+1
C      NtupleTag(m) = 'Pg12'
C      m = m+1
C      NtupleTag(m) = 'Pg13'
C      m = m+1
C      NtupleTag(m) = 'Pg14'
C      m = m+1
C      NtupleTag(m) = 'Pg21'
C      m = m+1
C      NtupleTag(m) = 'Pg22'
C      m = m+1
C      NtupleTag(m) = 'Pg23'
C      m = m+1
C      NtupleTag(m) = 'Pg24'
C      m = m+1
C      NtupleTag(m) = 'Pe11'
C      m = m+1
C      NtupleTag(m) = 'Pe12'
C      m = m+1
C      NtupleTag(m) = 'Pe13'
C      m = m+1
C      NtupleTag(m) = 'Pe14'
      m = m+1
      NtupleTag(m) = 'cutg'
      m = m+1
      NtupleTag(m) = 'cute'
      m = m+1
      NtupleTag(m) = 'cutr'
C      m = m+1
C      NtupleTag(m) = 'e1cut'
C      m = m+1
C      NtupleTag(m) = 'e2cut'
C      m = m+1
C      NtupleTag(m) = 'p1cut'
C
*
      NtupleSize = m
*
      call HBOOKN(id,title,NtupleSize,name,bank,NtupleTag) !create Ntuple
*
      call HCDIR(NtupleDirectory,'R') !record Ntuple directory
*
      call HCDIR(directory,' ')       !reset CERNLIB directory
*
      return
      end


