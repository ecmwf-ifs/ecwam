      SUBROUTINE MAKEGRID (BLOCK, GRID, IJS, IJL, IG)

! ----------------------------------------------------------------------

!**** *MAKEGRID* - MAKE GRIDDED WAM MODEL FIELDS FROM BLOCKED FIELDS.

!     H. GUNTHER       ECMWF    NOVEMBER 1989

!*    PURPOSE.
!     --------

!       GRIDDED WAVE FIELDS ARE CREATED FROM BLOCKED WAVE FIELDS.

!**   INTERFACE.
!     ----------

!        *CALL MAKEGRID*
!           *BLOCK*   REAL   DATA IN BLOCKED FORMAT
!           *GRID*    REAL   DATA IN GRID FORMAT
!           *IJS*     INTEGER INDEX OF FIRST GRIDPOINT IN BLOCK
!           *IJL*     INTEGER INDEX OF LAST  GRIDPOINT IN BLOCK
!           *IG*      INTEGER BLOCK NUMBER

!     METHOD.
!     -------

!       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
!       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
!       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
!       INITIALISED WITH ZMISS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWMAP   , ONLY : IXLG     ,KXLT




! ----------------------------------------------------------------------

      DIMENSION BLOCK(NIBLO), GRID(NGX,NGY)

! ----------------------------------------------------------------------

!*    1. INITIALIZATION
!     -----------------

      IF (IG.EQ.1) THEN
        DO J = 1,NGY
          DO I = 1,NGX
            GRID(I,J) = ZMISS
          ENDDO
        ENDDO
      ENDIF

!*    2. MAKE GRIDDED FIELD
!     ---------------------

      DO IJ = IJS,IJL
        IX = IXLG(IJ,IG)
        IY = NGY- KXLT(IJ,IG) +1
        GRID(IX,IY) = BLOCK(IJ)
      ENDDO

      RETURN
      END SUBROUTINE MAKEGRID
