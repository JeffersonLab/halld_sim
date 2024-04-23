      REAL FUNCTION GET_KIN(IFL,ITR)
C
C ===  Gets the parameters of a particle ITR
C===    IFL=0 - full momentum
C===       =1 - energy
C===       =2 - cos(THETA)
C===       =3 - TAN(theta)
C===       =4 - THETA (deg)
C===       =5 - PHI (deg)
C
      IMPLICIT NONE
      INTEGER  IFL,ITR
      INCLUDE ?
C
      INTEGER i
      REAL ptot,en,res,thet,sl,pt,phi,ct
C
C     -----------------------------------------------------------------
C
C
      GET_KIN=-9999.
      IF(IFL.LT.0.OR.IFL.GT.5) GO TO 999
      IF(ITR.LT.1.OR.ITR.GT.np) GO TO 999
C
      ptot=SQRT(pout(1,ITR)**2+pout(2,ITR)**2+pout(3,ITR)**2)
C
      en=SQRT(ptot**2+am(ITR)**2)
      ct=pout(3,ITR)/ptot
      pt=SQRT(pout(1,ITR)**2+pout(2,ITR)**2)
      sl=pt/pout(3,ITR)
      thet=ATAN2(pt,pout(3,ITR))
      phi=ATAN2(pout(2,ITR),pout(1,ITR))
C
      IF(IFL.EQ.0) THEN
         res=ptot
      ELSE IF(IFL.EQ.1) THEN
         res=en
      ELSE IF(IFL.EQ.2) THEN
         res=ct
      ELSE IF(IFL.EQ.3) THEN
         res=sl
      ELSE IF(IFL.EQ.4) THEN
         res=thet*180./3.1416
      ELSE IF(IFL.EQ.5) THEN
         res=phi*180./3.1416
         IF(res.LT.-90.) res=360.+res
      ENDIF
      GET_KIN=res
C
 999  RETURN
C
      END


