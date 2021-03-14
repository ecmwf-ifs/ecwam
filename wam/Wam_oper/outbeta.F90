      SUBROUTINE OUTBETA (IJS, IJL, PRCHAR, U10, US, Z0, CICVR, BETA)

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

!       *CALL* *OUTBETA (IJS, IJL, PRCHAR, U10, US, Z0, CICVR, BETA)
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *PRCHAR* - DEFAULT VALUE FOR CHARNOCK
!         *U10*    - WIND SPEED IN M/S.
!         *US*     - FRICTION VELOCITY IN M/S.
!         *Z0*     - ROUGHNESS LENGTH IN M.
!         *CICVR*  - SEA ICE COVER.
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

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOWPHYS  , ONLY : RNUM     ,ALPHAMIN  , ALPHAMAX
      USE YOWICE   , ONLY : LICERUN  ,LWAMRSETCI, CITHRSH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: PRCHAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10, US, Z0
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: BETA

      REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
      REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB
      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: Z0VIS, ZN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL)  :: ALPHAMAXU10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

!*    COMPUTE CHARNOCK 'CONSTANT' BETA.
!     ---------------------------------

      IF (LLGCBZ0) THEN
        ZN = RNUM
!!debile
        ZN = 0.0_JWRB
      ELSE
        ZN = 0.0_JWRB
      ENDIF
      ALPHAMAXU10(:)=MIN(ALPHAMAX,AMAX+BMAX*U10(:))

      DO IJ = IJS,IJL
        Z0VIS =  ZN/MAX(US(IJ),EPSUS)
        BETA(IJ) = G*(Z0(IJ)-Z0VIS)/MAX(US(IJ)**2,EPSUS)
        BETA(IJ) = MAX(MIN(BETA(IJ),ALPHAMAXU10(IJ)),ALPHAMIN)
      ENDDO

!      IF (LICERUN .AND. LWAMRSETCI) THEN
!        DO IJ = IJS,IJL
!          BETA(IJ) = (1.0_JWRB-CICVR(IJ))*BETA(IJ) + CICVR(IJ)*PRCHAR
!        ENDDO
!      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)

! ----------------------------------------------------------------------
      END SUBROUTINE OUTBETA
