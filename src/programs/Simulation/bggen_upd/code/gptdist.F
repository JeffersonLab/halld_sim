      REAL FUNCTION GPTDIST(TT,PP,NPAR)
C
C---      Returns the t-distribution - the value of a function depending on the coefficients PP
C
      IMPLICIT NONE
      REAL TT       ! (input) value of the argument
     +    ,PP(*)    ! (input) an array of coefficients, PP(k) 0<k<=NPAR
      INTEGER NPAR  ! max number of coefficients
C
C--- Let a1=PP(1), a2=PP(2) ...
C---  Function = a1 * (exp(t*a2)+exp(t*a3)*a4) * 1/(1-t*a5)**4
C--- for a4=a5=0 a single exponetial is left. 
C--- The normalization factor a1 is arbitrary (positive)
C
      REAL ff
C
      ff=PP(1)                                                     ! a scale
      IF(NPAR.GE.4) ff=ff*(EXP(TT*PP(2)) + EXP(TT*PP(3))*PP(4))   ! sum of 2 exponentials 
      IF(NPAR.GE.5) ff=ff*1./(1-TT*PP(5))**4                      ! dipole formula (from Frankfurt and Strikman)
C
      GPTDIST=ff
C
      RETURN
      END

