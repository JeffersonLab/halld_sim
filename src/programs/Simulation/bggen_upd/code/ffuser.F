      SUBROUTINE FFUSER(CKEY)
C
C
C---      Initialize the FFREAD and the relevant variables
C
      IMPLICIT NONE
      CHARACTER CKEY*16
C
      INCLUDE 'bg_partc.inc'
      INCLUDE 'bg_tmp.inc'
      INCLUDE 'bg_ctrl.inc'
C
      INTEGER kf,ig,i,j,ig1,lent
      CHARACTER cnam*20
C
C      write(6,*) ' FFUSER: ',CKEY,(K_PTMP(i),i=1,3)
      IF(CKEY.EQ.'PARTICLE') THEN
         ig=K_PTMP(1)
         kf=K_PTMP(2)
         IF(ig.GT.0.AND.ig.LE.MXPART) THEN
            IF(KF_PART(ig).EQ.0) THEN               
               NPADEF=NPADEF+1
               IDGDEF(NPADEF)=ig
            ENDIF
            IF(kf.EQ.0) kf=99999 ! not defined - should be non zero 
            DO j=1,NPADEF        ! Check if this KF is already used (apart from the current ig)
               ig1=IDGDEF(j)
               IF(ig1.NE.ig) THEN
                  IF(kf.EQ.KF_PART(ig1)) THEN
                     WRITE(6,1030) CKEY,(K_PTMP(i),i=1,3),ig1
 1030                FORMAT(' *** FFUSER error : ',A,3X,3I7
     +                ,' KF already set for IDG=',I4)
                  ENDIF
               ENDIF
            ENDDO
            KF_PART(ig)=kf
            LD_PART(ig)=K_PTMP(3)
            AM_PART(ig)=R_PTMP(4)
            WG_PART(ig)=R_PTMP(5)
            DO i=1,4
               KD_PART(i,ig)=K_PTMP(5+i)
            ENDDO
            CALL UHTOC(K_PTMP(10),4,cnam,16)
            lent=LEN_TRIM(cnam)
            CH_PART(ig)=cnam(1:lent)
         ELSE
            WRITE(6,1050) CKEY,(K_PTMP(i),i=1,3)
 1050       FORMAT(' *** FFUSER error: ',A,3X,3I7)
         ENDIF 
         DO i=1,20              !  reset the temp. arrays
            K_PTMP(i)=0
         ENDDO
         DO i=1,4
            CALL UCTOH('    ',K_PTMP(9+i),4,4)
         ENDDO
      ELSE IF(CKEY.EQ.'PYTMOD') THEN     ! PYTHIA input data, say 'PARP(2)=3.'
         IF(NPYMOD.LT.MXPYMOD) THEN
            NPYMOD=NPYMOD+1
            CALL UHTOC(KPYTMP,4,CPYMOD(NPYMOD),32)
         ELSE
            WRITE(6,1100) (KPYTMP(i),i=1,8), NPYMOD,MXPYMOD
 1100       FORMAT(' *** ERROR in ffuser.F: card *PYTMOD ',8A4
     +            ,' is ignored - too many of them, the limit '
     +            ,'MXPYMOD=',I5)
         ENDIF
         DO i=1,8
            CALL UCTOH('    ',KPYTMP(i),4,4)
         ENDDO         
      ENDIF
C
C
      RETURN
      END
