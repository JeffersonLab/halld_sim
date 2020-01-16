*
      double precision  a_0, r_p, A_n, sigma_by2, alpha
      common /nucl_par/ a_0, r_p, A_n, sigma_by2, alpha
      save   /nucl_par/
      double precision  rp_a0, ra_0, ro_0a, c
      common /nucl_nrm/ rp_a0, ra_0, ro_0a, c
      save   /nucl_nrm/
*
      double precision  ro_0
      common /obs_par/  ro_0
      save   /obs_par/
*
      double precision PI, DEGRAD
      parameter(PI = 3.141592653589793238d0)
      parameter(DEGRAD = PI / 180.d0)
*
      double precision PI05
      parameter(PI05 = 1.7724538509055159d0)      ! sqrt(PI)
*
