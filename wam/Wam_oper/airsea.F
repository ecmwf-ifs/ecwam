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

!       CALL TAUT_Z0

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : XKAPPA   ,XNLEV
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWWIND  , ONLY : WSPMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, KLEV, ICODE_WND

      REAL,DIMENSION(IJS:IJL),INTENT(IN)  ::  TAUW
      REAL,DIMENSION(IJS:IJL),INTENT(INOUT) :: U10, US
      REAL,DIMENSION(IJS:IJL),INTENT(OUT) ::  Z0

      INTEGER :: IJ, I, J

      REAL :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1 
      REAL :: XKAPPAD, XLOGLEV
      REAL :: XLEV
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)
#endif

!*    2. DETERMINE TOTAL STRESS (if needed)
!        ----------------------------------

      IF(ICODE_WND == 3 ) THEN
         XLEV = XNLEV(KLEV)
         CALL TAUT_Z0(IJS, IJL, XLEV, U10, TAUW, US, Z0)

      ELSE IF(ICODE_WND == 1 .OR. ICODE_WND == 2 ) THEN

!*    3. DETERMINE ROUGHNESS LENGTH (if needed).
!        ---------------------------

        CALL Z0WAVE (IJS, IJL, US, TAUW, Z0)

!*    3. DETERMINE U10 (if needed).
!        ---------------------------

        XKAPPAD=1.0/XKAPPA
        XLOGLEV=LOG(XNLEV(KLEV))

        DO IJ=IJS,IJL
          U10(IJ) = XKAPPAD*US(IJ)*(XLOGLEV-LOG(Z0(IJ)))
          U10(IJ) = MAX(U10(IJ),WSPMIN)
        ENDDO

      ELSE
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) ' + AIRSEA : INVALID VALUE OF ICODE_WND    +'
        WRITE(IU06,*) ' ICODE_WND = ',ICODE_WND
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE AIRSEA
