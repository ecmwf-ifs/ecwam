      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL, KLEV)

! ----------------------------------------------------------------------

!**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
!                                                   SQRT(TAUT)

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL STRESS.

!**   INTERFACE.
!     ----------

!       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL)*
!          *U10*  - INPUT BLOCK OF WINDSPEEDS U10.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *KLEV* - LEVEL HEIGHT INDEX

!     METHOD.
!     -------

!       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : XKAPPA   ,XNLEV
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,TAUT     ,DELTAUW  ,DELU

! ----------------------------------------------------------------------

      REAL,DIMENSION(IJS:IJL) ::  U10,TAUW,US,Z0

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

      ILEV=KLEV


!*    2. DETERMINE TOTAL STRESS FROM TABLE.
!        ----------------------------------

      DO IJ=IJS,IJL
        XI      = TAUW(IJ)/DELTAUW
        I       = MIN ( ITAUMAX-1, INT(XI) )
        DELI1   = MIN(1.,XI - REAL(I))
        DELI2   = 1. - DELI1
        XJ      = U10(IJ)/DELU
        J       = MIN ( JUMAX-1, INT(XJ) )
        DELJ1   = MIN(1.,XJ - REAL(J))
        DELJ2   = 1. - DELJ1
        US(IJ)  = (TAUT(I,J,ILEV)*DELI2 + TAUT(I+1,J,ILEV)*DELI1)*DELJ2
     &       +(TAUT(I,J+1,ILEV)*DELI2 + TAUT(I+1,J+1,ILEV)*DELI1)*DELJ1
      ENDDO

!*    3. DETERMINE ROUGHNESS LENGTH.
!        ---------------------------

      DO IJ=IJS,IJL
        SQRTCD  = MIN(U10(IJ)/US(IJ),100.0)
        Z0(IJ)  = XNLEV(ILEV)*EXP(-XKAPPA*SQRTCD)
      ENDDO

      RETURN
      END SUBROUTINE AIRSEA
