      REAL FUNCTION EFM_TYP(IG1,IG2)
C
C---      Eff mass of the particles of GEANT type IG1 and IG2
C
      IMPLICIT NONE
      INTEGER IG1,IG2
C
      INCLUDE ?
C
      INTEGER nm,im(4),i,j,k,ip
      REAL bm(4),pm(4,4),pouta(4),ef
C
      EFM_TYP=0.
      nm=0
      DO ip=1,NP
C         write(6,*) ip,ITYP(1,ip),IG1
         IF(ITYP(1,ip).EQ.IG1) THEN
            nm=nm+1
            bm(nm)=AM(ip)
            DO j=1,4
               pm(j,nm)=POUT(j,ip)
            ENDDO
         ENDIF
      ENDDO
      DO ip=1,NP
         IF(ITYP(1,ip).EQ.IG2) THEN
            nm=nm+1
            bm(nm)=AM(ip)
            DO j=1,4
               pm(j,nm)=POUT(j,ip)
            ENDDO
         ENDIF
      ENDDO
C
      IF(nm.EQ.2) THEN
         CALL EFMASS(nm,bm(1),pm(1,1),ef,pouta)
         EFM_TYP=ef
      ENDIF
C
      END
C
      INCLUDE 'efmass.f'
