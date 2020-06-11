C hbook.inc
C Sets up common block for hbook.
*
*-sizes of CERNLIB working space
*
      INTEGER iHbookSize,iHigzSize,iKuipSize,iPawSize
      PARAMETER (iHbookSize = 100000)
      Parameter (iHigzSize =  50000)
      PARAMETER (iKuipSize =  75000)
      PARAMETER (iPawSize = 3250000)
*
*-CERNLIB working space
*
        integer iCernMemory(iPawSize)
        COMMON /PAWC/ iCernMemory      !special nonstandard name!
*
C Ntuple ID stuff
*
	integer*4 		defaultID
	parameter     		(defaultID = 21)
      	character*132 		NtupleDirectory
        character*16            NtupleName
        integer*4               NtupleID,NtupleIO
  
        common                  NtupleID, NtupleDirectory, NtupleName, NtupleIO



