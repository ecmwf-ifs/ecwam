      SUBROUTINE MAKEGRID (BLOCK, GRID, IJS, IJL, IG, LSQRT, PMISS)

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
!           *LSQRT*   LOGICAL IF TRUE THE 4*SQUARE ROOT OF BLOCK IS
!                     COMPUTED WHEN TRANSFERRING TO GRID EXCEPT IF
!                     BLOCK EQUALS PMISS.
!           *PMISS*   REAL  VALUE GIVEN FOR ALL NON SEA POINTS. 

!     METHOD.
!     -------

!       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
!       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
!       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
!       INITIALISED WITH PMISS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO

! ----------------------------------------------------------------------

      LOGICAL LSQRT
      DIMENSION BLOCK(NIBLO), GRID(NGX,NGY)

!*    1. INITIALIZATION
!     -----------------

      CALL GSTATS(1996,0)
      IF (IG.EQ.1) THEN
        DO J = 1,NGY
          DO I = 1,NGX
            GRID(I,J) = PMISS
          ENDDO
        ENDDO
      ENDIF

!*    2. MAKE GRIDDED FIELD
!     ---------------------

      IF(LSQRT) THEN
        DO IJ = IJS,IJL
          IF(BLOCK(IJ).NE.PMISS) THEN
            IX = IXLG(IJ,IG)
            IY = NGY- KXLT(IJ,IG) +1
            GRID(IX,IY) = 4.*SQRT(BLOCK(IJ))
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
          GRID(IX,IY) = BLOCK(IJ)
        ENDDO
      ENDIF
      CALL GSTATS(1996,1)

      RETURN
      END SUBROUTINE MAKEGRID
