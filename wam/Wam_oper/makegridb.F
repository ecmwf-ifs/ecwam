      SUBROUTINE MAKEGRIDB (BLOCK, GRID, IJS, IJL, IG,
     &                     NGX, NGY, NIBLO, NBLO,
     &                     IXLG, KXLT)

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
!       INITIALISED WITH ZERO.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      PARAMETER(ZMISS=-999.)

! ----------------------------------------------------------------------

      DIMENSION BLOCK(NIBLO), GRID(NGX,NGY)
      DIMENSION IXLG(NIBLO,NBLO), KXLT(NIBLO,NBLO)
      LOGICAL LICERUN ; DATA LICERUN / .FALSE. /

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

!*    3. SET ICE POINTS TO MISSING DATA.
!     ----------------------------------

!     IF (LICERUN) THEN
!        DO 3100 KICE = 1,IGNO
!           IF (IGICE(KICE).EQ.IG) THEN
!              DO 3000 I = KSTART(KICE),KEND(KICE)
!                 IX = IXLG(ICEIJ(I),IG)
!                 IY = NGY- KXLT(ICEIJ(I),IG) +1
!                 GRID(IX,IY) = ZMISS
!           ENDIF
!     ENDIF

      RETURN
      END SUBROUTINE MAKEGRIDB




