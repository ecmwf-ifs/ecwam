      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL, KLEV)

! ----------------------------------------------------------------------

!**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
!                                                   SQRT(TAUT)
!     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
!                                                 TAUW

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

      USE YOWCOUP  , ONLY : ALPHA    ,XKAPPA   ,XNLEV
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,TAUT     ,DELTAUW     ,
     &              DELU   ,EPS1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, KLEV
      INTEGER :: ILEV, IJ, I, J
      REAL,DIMENSION(IJS:IJL),INTENT(IN)  ::  U10,TAUW
      REAL,DIMENSION(IJS:IJL),INTENT(OUT) ::  US,Z0
      REAL :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1 
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)
#endif

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

      ILEV=KLEV

!*    2. DETERMINE TOTAL STRESS FROM TABLE.
!        ----------------------------------

      DO IJ=IJS,IJL
        XI      = SQRT(TAUW(IJ))/DELTAUW
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

      CALL Z0WAVE (IJS, IJL, US, TAUW, Z0)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE AIRSEA
