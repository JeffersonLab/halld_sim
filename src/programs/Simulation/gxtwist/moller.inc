      real motar_angle_deg, motar_thick_cm, motar_z_cm
      real motar_zlimit_cm(2), motar_xydim_cm(2)
      integer motar_imate, motar_itmed, motar_ivolu, motar_irot
      common /moller_target/ motar_angle_deg,
     +                       motar_thick_cm,
     +                       motar_z_cm,
     +                       motar_zlimit_cm,
     +                       motar_xydim_cm,
     +                       motar_imate,
     +                       motar_itmed,
     +                       motar_ivolu,
     +                       motar_irot
      real mollerpar(3)
      equivalence (motar_angle_deg, mollerpar(1))

