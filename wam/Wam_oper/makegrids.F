      SUBROUTINE makegrids (kuso, pblock, psec4, kdl, kjs, kjl, kg,
     &                      klonrgg, kgx, kgy, kang, kfre, kxlg, kxlt,
     &                      kblo, kiblo)

!======================================================================

!**** *MAKEGRID* - MAKE GRIDDED WAM MODEL SPECTRA FROM BLOCKED SPECTRA.

!     B. HANSEN        ECMWF    MAY      1993

!*    PURPOSE.
!     --------

!       GRIDDED WAVE SPECTRA ARE CREATED FROM BLOCKED WAVE SPECTRA.

!**   INTERFACE.
!     ----------

!        *CALL MAKEGRID*
!           *KUSO*    REAL   DATA IN BLOCKED FORMAT
!           *PBLOCK*  REAL   DATA IN BLOCKED FORMAT
!           *PSEC4*   REAL   DATA IN GRID FORMAT
!           *KDL*     SIZE OF PSEC4
!           *KJS*     INTEGER INDEX OF FIRST GRIDPOINT IN BLOCK
!           *KJL*     INTEGER INDEX OF LAST  GRIDPOINT IN BLOCK
!           *KG*      INTEGER BLOCK NUMBER
!           *KLONRGG* INTEGER BLOCK NUMBER
!           *KGX*     TOTAL LENGTH OF FIELD IN X DIRECTION.
!           *KGY*     TOTAL LENGTH OF FIELD IN Y DIRECTION.
!           *KANG*    NUMBER OF DIRECTIONS.
!           *KFRE*    NUMBER OF FREQUENCIES.
!           *KXLG*    LONG. GRID INDEX.
!           *KXLT*    LAT. GRID INDEX.
!           *KBLO*    NUMBER OF BLOCKS.
!           *KIBLO*   NUMBER OF POINTS IN BLOCK.

!     METHOD.
!     -------

!       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
!       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
!       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
!       INITIALISED WITH ZERO.
!       IN THIS VERSION THE FIRST VALUE IN THE ARRAY IS LOCATED AT THE
!       NORTHENMOST AND WESTERNMOST CORNER OF THE AREA. THE NEXT POINT
!       LIES EAST OF THE FIRST POINT AND SO ON UNTIL THE END OF THE 
!       LATITUDE. THE THEN FOLLOWING VALUE LIES SOUTH OF THE FIRST
!       VALUE AND SO FORTH.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      INTEGER kxlg(kiblo,kblo), kxlt(kiblo,kblo), klonrgg(kgy)
      REAL    pblock (0:kiblo,kang,kfre)
      REAL    psec4 (kdl)
      INTEGER i_off (kjs:kjl)

! ----------------------------------------------------------------------

!* 0.    COMMON OUTPUT CONTROL.
!  ----------------------------

#include <outctrl.h>
      DATA csubna /"MAKEGRI"/

! ----------------------------------------------------------------------

!*    1. INITIALIZATION
!     -----------------

      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' SUB: ',csubna,': '

      IF (llsuvi) THEN
        WRITE(kuso,*) ' DO MAKEGRIDS FOR BLOCK ', kg
        WRITE(kuso,*) ' FIRST INDEX IN BLOCK   ', kjs
        WRITE(kuso,*) ' LAST INDEX IN BLOCK    ', kjl
      ENDIF

      ipsectrum = kang*kfre

      IF (lltest) WRITE(kuso,'(" ij  ix  iy klonrgg   io   i_off")')
!                               ....;...;...;......;......;......
      DO ij = kjs, kjl
        ix = kxlg(ij,kg)
        iy = kgy - kxlt(ij,kg)

        io = 0
!        DO ji = 1, iy
        DO ji = kgy, kxlt(ij,kg)+1, -1
          io = io + klonrgg(ji)
        ENDDO
        i_off(ij) = (io + ix-1) * ipsectrum

        IF (lltest .AND. iy .GT. 0)
     &   WRITE(kuso,'(3i4, 3i7)')ij, ix, iy, klonrgg(iy), io, i_off(ij)
        IF (lltest .AND. iy .LE. 0)
     &   WRITE(kuso,'(3i4, 3i7)')ij, ix, iy, 0,           io, i_off(ij)
      ENDDO

!*    2. MAKE GRIDDED FIELD
!     ---------------------

      DO jang = 1, kang
        DO jfre = 1, kfre
          jx0 = (jfre-1)*kang+1+(jang-1)
          DO ij = kjs, kjl
            jx = jx0 + i_off(ij)
            psec4(jx) = pblock(ij,jang,jfre)
          ENDDO
        ENDDO
      ENDDO

      msub = msub - 1
      RETURN
      END SUBROUTINE makegrids


