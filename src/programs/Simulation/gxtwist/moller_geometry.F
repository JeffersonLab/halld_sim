*
* moller_geometry - tweaks to the gxtwist tagger halld geometry
*                   for simulating the Hall D tagger used as a
*                   Moller spectrometer.
*
* author: richard.t.jones at uconn.edu
* version: february 10, 2024
*

      subroutine Moller_geometry_mods

#include "moller.inc"
#include "geant321/gcnum.inc"

      character*20 natmed
      character*4 chname,chshap,chmoth,chonly
      integer itmed,imat,isvol,ifield
      real fieldm,tmaxfd,stemax,deemax
      real epsil,stmin
      real nwbuf,ubuf(1)
      real boxpar(3)
      integer nboxpar
      integer ivolu
      integer irot
      real theta(3), phi(3)
      real xmotar,ymotar,zmotar

      itmed = NTMED + 1
      natmed = 'Moller target iron'
      imat = 10
      isvol = 0
      ifield = 1
      fieldm = 25
      tmaxfd = 1
      stemax = 1
      deemax = 0
      epsil = 0.001
      stmin = 0
      nwbuf = 0
      call gstmed(itmed,natmed,imat,isvol,ifield,fieldm,tmaxfd,
     +            stemax,deemax,epsil,stmin,ubuf,nwbuf)

      chname = 'FEMO'
      chshap = 'BOX '
      boxpar(1) = motar_xydim_cm(1) / 2.
      boxpar(2) = motar_xydim_cm(2) / 2.
      boxpar(3) = motar_thick_cm / 2.
      nboxpar = 3
      call gsvolu(chname,chshap,itmed,boxpar,nboxpar,ivolu)

      if (motar_z_cm .lt. motar_zlimit_cm(1) .or.
     +    motar_z_cm .gt. motar_zlimit_cm(2)) then
        print *, 'Error in control card FEMO: ',
     +           'third argument motar_z_cm =', motar_z_cm,
     +           'lies outside the prescribed bounds',
     +           motar_zlimit_cm(1), 'and', motar_zlimit_cm(2)
        stop
      endif

      irot = NROTM + 1
      theta(1) = 90
      phi(1) = 0
      theta(2) = 90 - motar_angle_deg
      phi(2) = 90
      theta(3) = motar_angle_deg
      phi(3) = 180
      call gsrotm(irot,theta(1),phi(1),theta(2),phi(2),theta(3),phi(3))

      nr = 1
      chmoth = 'QUAV'
      xmotar = 0
      ymotar = 0
      zmotar = motar_z_cm - (motar_zlimit_cm(1) + motar_zlimit_cm(2))/2.
      chonly = 'ONLY'
      call gspos(chname,nr,chmoth,xmotar,ymotar,zmotar,irot,chonly)

      motar_imate = imat
      motar_itmed = itmed
      motar_ivolu = ivolu
      motar_irot = irot
      end

      subroutine Moller_target_report
#include "moller.inc"
#include "geant321/gcbank.inc"
      character*20 chnama
      real A,Z,dens,radl,absl,ubuf(99)
      
      call GFMATE(motar_imate,chnama,A,Z,dens,radl,absl,ubuf,nwbuf)
      print *, 'Moller target report:'
      print *, ' * material ', chnama
      print *, '    A =', A
      print *, '    Z =', Z
      print *, '    dens =', dens
      print *, '    radl =', radl
      print *, '    absl =', absl
      end
