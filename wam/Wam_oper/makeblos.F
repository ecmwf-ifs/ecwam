      SUBROUTINE makeblos (pfl3, psec4, kdl, kjs, kjl, kblock,
     .                     klonrgg, kgx, kgy, kang, kfre, kxlg, kxlt,
     .                     kblo, kiblo)
C
C ======================================================================
C
C**** *MAKEBLOS* -  TRANSFER SPECTRA FROM GRID TO BLOCK.
C
C     BJORN HANSEN        ECMWF     MAY  1993
C
C     PURPOSE.
C     --------
C
C       TRANSFER SPECTRA FROM GRID TO BLOCK.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MAKEBLO (pfl3, GRID, kjs, kjl, kblock)*
C
C         *pfl3*    REAL       BLOCKED DATA.
C         *PSEC4*   REAL      GRIDDED DATA.
C         *KDL*     INTEGER   SIZE OF PSEC4.
C         *kjs*     INTEGER   INDEX OF FIRST POINT IN BLOCK.
C         *kjl*     INTEGER   INDEX OF LAST POINT IN BLOCK.
C         *kblock*  INTEGER   BLOCK NUMBER.
C
C
C ----------------------------------------------------------------------
C
C* 0.    COMMON OUTPUT CONTROL.
C  ----------------------------
C
#include <outctrl.h>
      DATA csubna /"MAKEBLO"/
C
C ----------------------------------------------------------------------
C
      INTEGER kxlg(kiblo,kblo)
      INTEGER kxlt(kiblo,kblo)
      INTEGER klonrgg(kgy)
      INTEGER i_off (kjs:kjl)
C
C ----------------------------------------------------------------------
C
      REAL pfl3(0:kiblo,kang,kfre)
      REAL psec4(kdl)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATIONS.
C        ----------------
 1000 CONTINUE
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(*,*)' ',csubna,': BLOCK=',kblock
C
      IF (lltest)
     .write(81,'(a40)' )'       JX,  KGX,  KGY, KFRE,  KANG,    JI'

      ipsectrum = kang*kfre

      DO ji = kjs, kjl
         ix = kxlg(ji,kblock)
         iy = kgy - kxlt(ji,kblock)

         io = 0
C        DO jy = 1, iy
         DO jy = kgy, kxlt(ji,kblock)+1, -1
           io = io + klonrgg(jy)
         ENDDO
         i_off(ji) = (io + ix-1) * ipsectrum

      ENDDO


C
C*    2. MAKE BLOCKED FIELD
C     ---------------------
C

      DO 2300 jang = 1, kang
        DO 2200 jfre = 1, kfre
          jx0 = (jfre-1)*kang+1+(jang-1)
          DO 2100 ji = kjs, kjl
            jx = jx0 + i_off(ji)
            pfl3(ji,jang,jfre) = psec4(jx)
            IF (lltest)
     .      WRITE (81,1300) jx, ix, iy, jfre, jang, ji, psec4(jx)
 2100     CONTINUE
 2200   CONTINUE
 2300 CONTINUE

C     DO 1200 ji = kjs,kjl
C        ix = kxlg(ji,kblock)
C        iy = kgy-kxlt(ji,kblock)
C        i_off = ispoff(klonrgg, kgy, kang, kfre, iy, ix)
C        DO jang=1,kang
C           DO jfre=1,kfre
C              jx=(jfre-1)*kang+1+(jang-1)+i_off
C              pfl3(ji,jang,jfre) = psec4(jx)
C              IF (lltest)
C    .         write(81,1300) jx, ix, iy, jfre, jang, ji, psec4(jx)
C           ENDDO
C        ENDDO
C1200 CONTINUE
 
      msub = msub - 1
      RETURN
 1300 FORMAT(i10,5i6,1x,f15.13)
      END
