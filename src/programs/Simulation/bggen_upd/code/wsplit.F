      SUBROUTINE WSPLIT(CLIN,NWMX,NW,IWOR,LWOR,ITYPW)
C
C---     Split a line CLIN into words CWOR
C
      IMPLICIT NONE
C--   INPUT:
      CHARACTER CLIN*(*) ! string
      INTEGER NWMX       ! max number of words
C--   OUTPUT:
      INTEGER NW         ! number of words
     +       ,IWOR(*)    ! the start address
     +       ,LWOR(*)    ! length of each word
     +       ,ITYPW(*)   ! type of each word: 0 - not a number, =1 integer, =2 - real
C
      CHARACTER cwor*160
C
      INTEGER il,il1,lenl,lw,i
C
      NW=0
      lenl=LEN(CLIN)
C
      il=1
C
 10   IF(CLIN(il:il).NE.' ') THEN
C
         IF(NW.LT.NWMX) THEN
            NW=NW+1
            IWOR(NW)=il
         ENDIF
         cwor=' '
         lw=INDEX(CLIN(il:lenl),' ')-1
         IF(lw.LE.0) lw=lenl-il+1
         cwor=CLIN(il:il+lw-1)
         il=il+lw
         LWOR(NW)=lw
         ITYPW(NW)=1
         DO i=1,lw
            IF(cwor(i:i).NE.'.'.AND.
     +         cwor(i:i).NE.'E'.AND.cwor(i:i).NE.'e'.AND.
     +       (LLT(cwor(i:i),'0').OR.LGT(cwor(i:i),'9'))) 
     +           ITYPW(NW)=0
         ENDDO
         IF(ITYPW(NW).EQ.1.AND.INDEX(cwor,'.').GT.0) ITYPW(NW)=2
C         WRITE(6,*) NW,LWOR(NW),' ',CWOR(NW)
      ELSE
         il=il+1
      ENDIF

      IF(il.LT.lenl) GO TO 10
C
      RETURN
      END
