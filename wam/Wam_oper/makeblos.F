      SUBROUTINE makeblos (pfl3, psec4, kdl, kjs, kjl, kblock,
     &                     klonrgg, kgx, kgy, kang, kfre, kxlg, kxlt,
     &                     kblo, kiblo)

! ======================================================================

!**** *MAKEBLOS* -  TRANSFER SPECTRA FROM GRID TO BLOCK.

!     BJORN HANSEN        ECMWF     MAY  1993

!     PURPOSE.
!     --------

!       TRANSFER SPECTRA FROM GRID TO BLOCK.

!**   INTERFACE.
!     ----------

!       *CALL* *MAKEBLO (pfl3, GRID, kjs, kjl, kblock)*

!         *pfl3*    REAL       BLOCKED DATA.
!         *PSEC4*   REAL      GRIDDED DATA.
!         *KDL*     INTEGER   SIZE OF PSEC4.
!         *kjs*     INTEGER   INDEX OF FIRST POINT IN BLOCK.
!         *kjl*     INTEGER   INDEX OF LAST POINT IN BLOCK.
!         *kblock*  INTEGER   BLOCK NUMBER.


! ----------------------------------------------------------------------

!* 0.    COMMON OUTPUT CONTROL.
!  ----------------------------

#include <outctrl.h>
      DATA csubna /"MAKEBLO"/

! ----------------------------------------------------------------------

      INTEGER kxlg(kiblo,kblo)
      INTEGER kxlt(kiblo,kblo)
      INTEGER klonrgg(kgy)
      INTEGER i_off (kjs:kjl)

! ----------------------------------------------------------------------

      REAL pfl3(0:kiblo,kang,kfre)
      REAL psec4(kdl)

! ----------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(*,*)' ',csubna,': BLOCK=',kblock

      IF (lltest)
     &write(81,'(a40)' )'       JX,  KGX,  KGY, KFRE,  KANG,    JI'

      ipsectrum = kang*kfre

      DO ji = kjs, kjl
        ix = kxlg(ji,kblock)
        iy = kgy - kxlt(ji,kblock)

        io = 0
!        DO jy = 1, iy
        DO jy = kgy, kxlt(ji,kblock)+1, -1
          io = io + klonrgg(jy)
        ENDDO
        i_off(ji) = (io + ix-1) * ipsectrum

      ENDDO



!*    2. MAKE BLOCKED FIELD
!     ---------------------


      DO jang = 1, kang
        DO jfre = 1, kfre
          jx0 = (jfre-1)*kang+1+(jang-1)
          DO ji = kjs, kjl
            jx = jx0 + i_off(ji)
            pfl3(ji,jang,jfre) = psec4(jx)
            IF (lltest)
     &      WRITE (81,1300) jx, ix, iy, jfre, jang, ji, psec4(jx)
          ENDDO
        ENDDO
      ENDDO

!     DO 1200 ji = kjs,kjl
!        ix = kxlg(ji,kblock)
!        iy = kgy-kxlt(ji,kblock)
!        i_off = ispoff(klonrgg, kgy, kang, kfre, iy, ix)
!        DO jang=1,kang
!           DO jfre=1,kfre
!              jx=(jfre-1)*kang+1+(jang-1)+i_off
!              pfl3(ji,jang,jfre) = psec4(jx)
!              IF (lltest)
!    .         write(81,1300) jx, ix, iy, jfre, jang, ji, psec4(jx)
!           ENDDO
!        ENDDO
!1200 CONTINUE
 
      msub = msub - 1
      RETURN
 1300 FORMAT(i10,5i6,1x,f15.13)
      END SUBROUTINE makeblos


