      SUBROUTINE makegrids (kuso, pblock, psec4, kdl, kjs, kjl, kg,
     .                      klonrgg, kgx, kgy, kang, kfre, kxlg, kxlt,
     .                      kblo, kiblo)
C
C======================================================================
C
C**** *MAKEGRID* - MAKE GRIDDED WAM MODEL SPECTRA FROM BLOCKED SPECTRA.
C
C     B. HANSEN        ECMWF    MAY      1993
C
C*    PURPOSE.
C     --------
C
C       GRIDDED WAVE SPECTRA ARE CREATED FROM BLOCKED WAVE SPECTRA.
C
C**   INTERFACE.
C     ----------
C
C        *CALL MAKEGRID*
C           *KUSO*    REAL   DATA IN BLOCKED FORMAT
C           *PBLOCK*  REAL   DATA IN BLOCKED FORMAT
C           *PSEC4*   REAL   DATA IN GRID FORMAT
C           *KDL*     SIZE OF PSEC4
C           *KJS*     INTEGER INDEX OF FIRST GRIDPOINT IN BLOCK
C           *KJL*     INTEGER INDEX OF LAST  GRIDPOINT IN BLOCK
C           *KG*      INTEGER BLOCK NUMBER
C           *KLONRGG* INTEGER BLOCK NUMBER
C           *KGX*     TOTAL LENGTH OF FIELD IN X DIRECTION.
C           *KGY*     TOTAL LENGTH OF FIELD IN Y DIRECTION.
C           *KANG*    NUMBER OF DIRECTIONS.
C           *KFRE*    NUMBER OF FREQUENCIES.
C           *KXLG*    LONG. GRID INDEX.
C           *KXLT*    LAT. GRID INDEX.
C           *KBLO*    NUMBER OF BLOCKS.
C           *KIBLO*   NUMBER OF POINTS IN BLOCK.
C
C     METHOD.
C     -------
C
C       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
C       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
C       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
C       INITIALISED WITH ZERO.
C       IN THIS VERSION THE FIRST VALUE IN THE ARRAY IS LOCATED AT THE
C       NORTHENMOST AND WESTERNMOST CORNER OF THE AREA. THE NEXT POINT
C       LIES EAST OF THE FIRST POINT AND SO ON UNTIL THE END OF THE 
C       LATITUDE. THE THEN FOLLOWING VALUE LIES SOUTH OF THE FIRST
C       VALUE AND SO FORTH.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      INTEGER kxlg(kiblo,kblo), kxlt(kiblo,kblo), klonrgg(kgy)
      REAL    pblock (0:kiblo,kang,kfre)
      REAL    psec4 (kdl)
      INTEGER i_off (kjs:kjl)
C
C ----------------------------------------------------------------------
C
C* 0.    COMMON OUTPUT CONTROL.
C  ----------------------------
C
#include <outctrl.h>
      DATA csubna /"MAKEGRI"/
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATION
C     -----------------
C
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' SUB: ',csubna,': '
C
      IF (llsuvi) THEN
        WRITE(kuso,*) ' DO MAKEGRIDS FOR BLOCK ', kg
        WRITE(kuso,*) ' FIRST INDEX IN BLOCK   ', kjs
        WRITE(kuso,*) ' LAST INDEX IN BLOCK    ', kjl
      ENDIF

      ipsectrum = kang*kfre

      IF (lltest) WRITE(kuso,'(" ij  ix  iy klonrgg   io   i_off")')
C                               ....;...;...;......;......;......
      DO ij = kjs, kjl
         ix = kxlg(ij,kg)
         iy = kgy - kxlt(ij,kg)

         io = 0
C        DO ji = 1, iy
         DO ji = kgy, kxlt(ij,kg)+1, -1
           io = io + klonrgg(ji)
         ENDDO
         i_off(ij) = (io + ix-1) * ipsectrum
C
         IF (lltest .AND. iy .GT. 0)
     &   WRITE(kuso,'(3i4, 3i7)')ij, ix, iy, klonrgg(iy), io, i_off(ij)
         IF (lltest .AND. iy .LE. 0)
     &   WRITE(kuso,'(3i4, 3i7)')ij, ix, iy, 0,           io, i_off(ij)
      ENDDO
C
C*    2. MAKE GRIDDED FIELD
C     ---------------------
C
      DO 2300 jang = 1, kang
        DO 2200 jfre = 1, kfre
          jx0 = (jfre-1)*kang+1+(jang-1)
          DO 2100 ij = kjs, kjl
            jx = jx0 + i_off(ij)
            psec4(jx) = pblock(ij,jang,jfre)
 2100     CONTINUE
 2200   CONTINUE
 2300 CONTINUE
C
      msub = msub - 1
      RETURN
      END
