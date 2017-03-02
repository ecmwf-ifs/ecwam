      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL, KLEV, ICODE_WND)

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

!       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL, KLEV, ICODE_WND)*
!          *U10*  - INPUT OR OUTPUT BLOCK OF WINDSPEED U10.
!          *TAUW* - INPUT BLOCK OF WAVE STRES.
!          *US*   - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *KLEV* - LEVEL HEIGHT INDEX
!          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN UPDATED:
!                     U10: ICODE_WND=3 --> US will be updated
!                     US:  ICODE_WND=1 OR 2 --> U10 will be updated

   

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
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,TAUT     ,DELTAUW     ,
     &              DELU   ,EPS1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWWIND  , ONLY : WSPMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, KLEV, ICODE_WND

      REAL,DIMENSION(IJS:IJL),INTENT(IN)  ::  TAUW
      REAL,DIMENSION(IJS:IJL),INTENT(INOUT) :: U10, US
      REAL,DIMENSION(IJS:IJL),INTENT(OUT) ::  Z0

      INTEGER :: ILEV, IJ, I, J

      REAL :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1 
      REAL :: XKAPPAD, XLOGLEV
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)
#endif

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

      ILEV=KLEV

!*    2. DETERMINE TOTAL STRESS FROM TABLE. (if needed)
!        ----------------------------------

      IF(ICODE_WND == 3 ) THEN
        DO IJ=IJS,IJL
          XI    = SQRT(TAUW(IJ))/DELTAUW
          I     = MIN ( ITAUMAX-1, INT(XI) )
          DELI1 = MIN(1.,XI - REAL(I))
          DELI2 = 1. - DELI1
          XJ    = U10(IJ)/DELU
          J     = MIN ( JUMAX-1, INT(XJ) )
          DELJ1 = MIN(1.,XJ - REAL(J))
          DELJ2 = 1. - DELJ1
          US(IJ)= (TAUT(I,J,ILEV)*DELI2 + TAUT(I+1,J,ILEV)*DELI1)*DELJ2
     &       +(TAUT(I,J+1,ILEV)*DELI2 + TAUT(I+1,J+1,ILEV)*DELI1)*DELJ1
        ENDDO
      ENDIF

!*    3. DETERMINE ROUGHNESS LENGTH.
!        ---------------------------

      CALL Z0WAVE (IJS, IJL, US, TAUW, Z0)


!*    3. DETERMINE U10 (if needed).
!        ---------------------------
      IF(ICODE_WND == 1 .OR. ICODE_WND == 2 ) THEN

        XKAPPAD=1.0/XKAPPA
        XLOGLEV=LOG(XNLEV(ILEV))

        DO IJ=IJS,IJL
          U10(IJ) = XKAPPAD*US(IJ)*(XLOGLEV-LOG(Z0(IJ)))
          U10(IJ) = MAX(U10(IJ),WSPMIN)
        ENDDO

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE AIRSEA
