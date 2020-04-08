      SUBROUTINE OUTBETA (IJS, IJL, U10, US, Z0, CICVR, PRCHAR, BETA)

! ----------------------------------------------------------------------

!**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING 
!     J.BIDLOT       ECMWF       MARCH 2014  USE MODEL CHARNOCK FOR ALL
!                                SEA POINTS (INCLUDING SEA ICE POINTS)

!*    PURPOSE.
!     --------

!       COMPUTES THE CHARNOCK PARAMETER.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (IJS, IJL, US, Z0, CICVR, BETA)
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *U10*    - WIND SPEED IN M/S.
!         *US*     - FRICTION VELOCITY IN M/S.
!         *Z0*     - ROUGHNESS LENGTH IN M.
!         *CICVR*  - SEA ICE COVER.
!         *PRCHAR* - DEFAULT VALUE FOR CHARNOCK
!         *BETA*   - CHARNOCK FIELD (BLOCK ARRAY)
!         

!     EXTERNALS.
!     ----------

!       NONE.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE  ,CITHRSH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10, US, Z0
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), INTENT(IN) :: PRCHAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: BETA

      REAL(KIND=JWRB), PARAMETER :: ALPHAMAX=0.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
      REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB
      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: ALPHAMAXU10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

      DO IJ = IJS,IJL
        BETA(IJ) = G*Z0(IJ)/MAX(US(IJ)**2,EPSUS)
        ALPHAMAXU10=MIN(ALPHAMAX,AMAX+BMAX*U10(IJ))
        BETA(IJ) = MIN(BETA(IJ),ALPHAMAXU10)
      ENDDO

      IF (LICERUN .AND. LMASKICE) THEN
        DO IJ = IJS,IJL
          BETA(IJ) = (1.0_JWRB-CICVR(IJ))*BETA(IJ) + CICVR(IJ)*PRCHAR
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)

! ----------------------------------------------------------------------
      END SUBROUTINE OUTBETA
