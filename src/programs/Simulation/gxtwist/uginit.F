*
* $Id$
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE UGINIT
*
************************************************************************
*                                                                      *
*              To initialise GEANT3 program and read data cards        *
*                                                                      *
************************************************************************
*
#include "geant321/gckine.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcphys.inc"
#include "hdtrackparams.inc"
#include "halo.inc"
#include "moller.inc"
*
      integer iskip
      integer infile(20)
      integer outfile(20)
      common /inputFile/ iskip,infile,outfile
      data iskip/0/
      data infile/20*0/
      data outfile/20*0/
      integer openInput, skipInput, openOutput
      external openInput, skipInput, openOutput
      real beamE0, beamEpeak, beamEmin, radColDist, colDiam
      real beamEmit, radThick
      common /beamPars/ beamE0,beamEpeak,beamEmin,radColDist,colDiam,
     +                  beamEmit, radThick
      data beamE0/0/
      data beamEpeak/0/
      data beamEmin/0/
      data radColDist/76.0/
      data colDiam/0.0034/
      data beamEmit/1e-8/
      data radThick/20e-6/

      data fhalo/0.0/
      data bfield/-999.0/
      data nosecondaries/0/

      data motar_angle_deg/0.0/
      data motar_thick_cm/0.0/
      data motar_z_cm/0.0/
      data motar_zlimit_cm/65.5,290.7/
      data motar_xydim_cm/1.0,1.0/

C  Use this parameter to set up a minimum photon energy 
C  for the coherent bremsstrahlung beam generator - see beamgen.F
      real    xMinimum,freqMaximum,beamStartZ,Theta02
      common /coherentGen/xMinimum,freqMaximum,beamStartZ,Theta02


*
*     -----------------------------------------------------------------
*
*             Initialize GEANT
C..geant..
      call trapfpe()
      CALL GINIT
*
*             Prints version number
*
      WRITE(LOUT,1000)
*
*             IKINE  = particle type (default=1=gamma)
*             PKINE(1)=particle energy
*             IKINE and PKINE can be changed with the data card KINE
*
      PKINE(1)=10.
      PKINE(5)=4.
      IKINE=1
      ICOMP=1
      IPAIR=1
      IBREM=1
      IANNI=1
      IPHOT=1
      IHADR=1
*
*             Initialize GEANT/ZBOOK data structures
*
C..geant..
      CALL GZINIT
*
*             Define user FFREAD data cards (format free input)
*
*
*             Read the data cards
*
      OPEN(UNIT=4,FILE='control.in',STATUS='UNKNOWN')
      CALL FFSET('LINP',4)
      call ffkey('infile',infile,20,'MIXED')
      call ffkey('skip',iskip,1,'INTEGER')
      call ffkey('outfile',outfile,20,'MIXED')
      call ffkey('beam',beamE0,7,'REAL')
      call ffkey('bfield', bfield,1,'REAL')
      call ffkey('nosecondaries', nosecondaries,1,'INTEGER')
      call ffkey('halo',fhalo,1,'REAL')
      call ffkey('femollertarget',mollerpar,3,'REAL')
      call gtgamaff()
      CALL GFFGO
*
*             Verify that the random number seeds are OK
*
      call GRNDMQ(iseed1,iseed2,0,'G')
      if ((iseed1.eq.0).and.(iseed2.eq.0)) then
        write(LOUT,980) 'UGINIT error: '
980     format(a14,'initial random number generator seeds are 0!')
        stop 'cannot continue without good random numbers'
      endif
*
*             Open the input stream
*
      if (infile(1) .gt. 0) then
        ifail = openInput(infile)
        if (ifail .ne. 0) then
           write(lout,9000) infile
 9000      format('GUKINE ERROR - Could not open input stream ',20a4)
           stop
        endif
        if (iskip .gt. 0) then
          ifail = skipInput(iskip)
        endif
        if (nevent .eq. 0) then
          nevent = 999999999
        endif
      endif
*
*             Open the output stream
*
      if (outfile(1) .eq. 0) then
        call UCTOH('gxtwist.hddm',outfile,4,12)
      endif
      ifail = openOutput(outfile)
      if (ifail .ne. 0) then
        write(lout,9010) outfile
 9010   format('GUKINE ERROR - Could not open output stream ',20a4)
        stop
      endif
*
*             Assign the beam parameters
*
      if (beamEmin.lt.0) then
         print *
         print *, 'Error in uginit:',
     +        ' beamEmin is specified with negative value,',
     +        ' cannot continue.'
         stop
      elseif (beamEmin.gt.0.12) then
         print *
         print *, 'Error in uginit:',
     +        ' beamEmin is larger than a default value of 0.12 GeV,',
     +        ' cannot continue.'
         stop
      elseif (beamEmin.lt.1e-6) then
         print *,
     +        ' beamEmin is smaller than a minimum value of 1 KeV,',
     +        ' Using default (120 MeV).'
         beamEmin = 0.12
      endif
      xMinimum = beamEmin/beamE0
            
      call cobrems(beamE0,beamEpeak,beamEmit,radThick,
     +             radColDist,colDiam,0)
*
*             Initialize graphics package
*
      CALL GDINIT
*
*             Initialize the Hall D geometry
*
      call HDDSgeant3
      if (motar_thick_cm .gt. 0) then
        call Moller_geometry_mods
      endif
      call Goptimize
*
*             Open the HBOOK file for output
*
      call HROPEN(50,'RZfile','geant.hbook','N',65536,istat)
*
*             Initialize GEANT tracking structures
*
C..geant..
      CALL GGCLOS
      CALL GPART
      CALL GPHYSI
      call Moller_target_report
      CALL GRFILE(3,'gxtwist.rz','ON')
      CALL HCDIR('//RZfile',' ')

*
*             Initialize the GELHAD package and verify parameters
*
      call gelh_vrfy()
*
*             Load FLUKA and MICAP cross section data if selected
*
      if (IHADR.ge.3) call FLINIT
      if (IHADR.eq.4) call GMORIN

 1000 FORMAT(/,'  MODE VERSION 1.00 : ',/)
      END
