      REAL FUNCTION GPTDIST(TTT) ! ,PP,NPAR)
C
C---      Returns the t-distribution - the value of a function depending on the coefficients PP
C
      REAL TT,TTT       ! (input) value of the argument
     +    ,PP(11)    ! (input) an array of coefficients, PP(k) 0<k<=NPAR
      INTEGER NPAR  ! max number of coefficients
C
C--- Let a1=PP(1), a2=PP(2) ...
C---  Function = a1 * (exp(t*a2)+exp(t*a3)*a4) * 1/(1-t*a5)**4
C--- for a4=a5=0 a single exponetial is left. 
C--- The normalization factor a1 is arbitrary (positive)
C
      REAL ff
C
C REACTDIS1  9.28  1.    1.52  -0.650 0.00119 
C REACTDIS2 10.36  1.    1.22   0.268 0.033   
C REACTDIS3 13.    1.    1.89   0.646 0.149   

      NPAR=10
      PP(1)=1.
      PP(2)=1.52
      PP(3)=-0.65
      PP(4)=0.00119
      PP(5)=0.

      PP(2)=1.22
      PP(3)=0.268
      PP(4)=0.033

      PP(2)=1.89
      PP(3)=0.646
      PP(4)=0.149

      PP(2)=1.52
      PP(3)=-0.65
      PP(4)=0.00119

      TT=-TTT
      ff=PP(1)                                                     ! a scale
      IF(NPAR.GE.4) ff=ff*(EXP(TT*PP(2)) + EXP(TT*PP(3))*PP(4))   ! sum of 2 exponentials 
      IF(NPAR.GE.5) ff=ff*1./(1-TT*PP(5))**4                      ! dipole formula (from Frankfurt and Strikman)
C
      GPTDIST=ff
C
      RETURN
      END

