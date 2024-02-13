C=======================================================================
      REAL FUNCTION SF_LOOKUP(PP,EE)
C
C----          Returns the spectral probability (SF*4pi*p**2*dP*dE)
C----          for momentum PP, energy EE
C
      IMPLICIT NONE
      REAL PP,EE
C
      INCLUDE 'sf_data.inc'
C
      INTEGER ip1,ie1,ip2,ie2
      REAL p1,p2,e1,e2,dp,de,wp,we,sf1,sf2,fun,sf(2,2)
C
C     ------------------------------------------------------------------
C
      SF_LOOKUP=0.
C
C---    Find the P cell
C
      IF(PP.LT.SF_PCEN(1)) THEN
         ip1=1
      ELSE IF(PP.GE.SF_PCEN(NP_SF)) THEN
         ip1=NP_SF-1
      ELSE
         ip1=2
         DO WHILE(ip1.LE.NP_SF.AND.PP.GE.SF_PCEN(ip1))
            ip1=ip1+1
         ENDDO
         ip1=ip1-1
      ENDIF
      ip2=ip1+1
      p1=SF_PCEN(ip1)
      p2=SF_PCEN(ip2)
      dp=(p2-p1)
      wp=PP-p1
C
C---    Find the E cell
C
      IF(EE.LT.SF_ECEN(1)) THEN
         ie1=1
      ELSE IF(EE.GE.SF_ECEN(NE_SF)) THEN
         ie1=NE_SF-1
      ELSE
         ie1=2
         DO WHILE(ie1.LE.NE_SF.AND.EE.GE.SF_ECEN(ie1))
            ie1=ie1+1
         ENDDO
         ie1=ie1-1
      ENDIF
      ie2=ie1+1
      e1=SF_ECEN(ie1)
      e2=SF_ECEN(ie2)
      de=e2-e1
      we=EE-e1
C
      sf(1,1)=SF_PROB(ie1,ip1)
      sf(1,2)=SF_PROB(ie1,ip2)
      sf(2,1)=SF_PROB(ie2,ip1)
      sf(2,2)=SF_PROB(ie2,ip2)
C
C---    Correction for EE<SF_ECEN(1) and sf rising toward e->0 - add sf=0 at e1=ecen(1)-dde
C---                             dde=MAX(ebin(i),5.) 
C
      IF(EE.LT.SF_ECEN(1)) THEN
         IF(SF_PROB(ie1,ip1).GT.SF_PROB(ie2,ip1)) THEN
            e2=e1
            e1=e2-MAX(SF_EWID(1)/2.,5.)
            IF(e1.LT.0.) e1=0.
            de=e2-e1
            we=EE-e1
            sf(1,1)=0.
            sf(1,2)=0.
            sf(2,1)=SF_PROB(1,ip1)
            sf(2,2)=SF_PROB(1,ip2)
C            write(6,*) EE,PP,ip1,e1,e2,SF_PROB(ie1,ip1),SF_PROB(ie2,ip1)
C     +                , sf(1,1),sf(2,1)
         ENDIF
      ENDIF
C
C---    Linear interpolation
C
      sf1=sf(1,1)+(sf(1,2)-sf(1,1))*wp/dp
      sf2=sf(2,1)+(sf(2,2)-sf(2,1))*wp/dp
      IF(sf1.LT.0.) sf1=0.
      IF(sf2.LT.0.) sf2=0.
C
      fun=sf1+(sf2-sf1)*we/de
C      IF(EE.LT.20..AND.PP.LT.50.) write(6,3000) PP,EE,p1,e1
C     +             ,ip1,ip2,ie1,ie2
C     +             ,wp/dp,we/de,sf,sf1,sf2,fun
C 3000 FORMAT(2F6.1,2X,2F6.2,1X,4I4,1X,2E11.2,1X,4E11.3,1X,3E13.4)
C
      IF(fun.LE.0.) GO TO 999
      SF_LOOKUP=fun
C
 999  CONTINUE
      RETURN
      END
