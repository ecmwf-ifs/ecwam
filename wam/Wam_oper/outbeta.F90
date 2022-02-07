SUBROUTINE OUTBETA (KIJS, KIJL,           & 
 &                  U10, USTAR, Z0M, Z0B, &
 &                  CHRNCK , BETAHQ, CD)

! ----------------------------------------------------------------------

!**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING 
!     J.BIDLOT       ECMWF       MARCH 2014  USE MODEL CHARNOCK FOR ALL
!                                SEA POINTS

!*    PURPOSE.
!     --------

!       COMPUTES THE CHARNOCK PARAMETER.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (KIJS, KIJL,
!                        U10, USTAR, Z0M, Z0B,
!                        CHRNCK , BETAHQ, CD)
!         *KIJS*    - INDEX OF FIRST GRIDPOINT.
!         *KIJL*    - INDEX OF LAST GRIDPOINT.
!         *U10*     - WIND SPEED.
!         *USTAR*   - FRICTION VELOCITY.
!         *Z0M*     - ROUGHNESS LENGTH.
!         *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
!         *BETA*    - CHARNOCK
!         *BETAHQ*  - EQUIVALENT CHARNOCK FIELD FOR HEAT AND MOISTURE
!                    (i.e. Charnock with the background roughness removed)   
!         *CD*      - OPTIONAL : DRAG COEFFICIENT THAT IS CONSISTENT WITH THE DERIVED CHARNOCK
!                     AND WHETHER OR NOT VISCOSITY EFFECT WAS ACCOUNTED FOR IN Z0M
!                     ASSUMING A LOGARTHMIC PROFILE FOR THE WIND.
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
      USE YOWPCONS , ONLY : G , GM1, EPSUS
      USE YOWPHYS  , ONLY : XKAPPA, XNLEV, RNUM , PRCHAR, ALPHAMIN, ALPHAMAX, ALPHA

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: U10
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: Z0B
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: CHRNCK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: BETAHQ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT), OPTIONAL :: CD


      REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
      REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB

!     !!!!! BETAHQ_REDUCE was introduced to reduce the impact of the sea state dependent
!     !!!!! heat and moisture surface flux. It is just a dirty tuning knob !!

      REAL(KIND=JWRB), PARAMETER :: BETAHQ_REDUCE=0.25_JWRB

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: ZN, GUSM2, Z0ATM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: USM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPHAMAXU10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)


!*    COMPUTE CHARNOCK 'CONSTANT' CHRNCK.
!     ---------------------------------

      IF (LLGCBZ0) THEN
        ZN = RNUM
        ALPHAMAXU10(:)=ALPHAMAX
      ELSE
        ZN = 0.0_JWRB
        ALPHAMAXU10(:)=MIN(ALPHAMAX,AMAX+BMAX*U10(:))
      ENDIF


      DO IJ = KIJS,KIJL
        USM(IJ) = 1.0_JWRB/MAX(USTAR(IJ), EPSUS)
        GUSM2 = G*USM(IJ)**2
        CHRNCK(IJ) = MAX(MIN(CHRNCK(IJ), ALPHAMAXU10(IJ)), ALPHAMIN)
        BETAHQ(IJ) = MAX( BETAHQ_REDUCE*(CHRNCK(IJ)-MAX(Z0B(IJ)*GUSM2, ALPHA)), ALPHAMIN)
      ENDDO

      IF( PRESENT(CD) ) THEN
        DO IJ = KIJS,KIJL
!!!     we are assuming here that z0 ~ ZN/USTAR + Charnock USTAR**2/g
!!!     in order to fit with what is used in the IFS.
          Z0ATM = ZN*USM(IJ) + GM1 * CHRNCK(IJ) * MAX(USTAR(IJ), EPSUS)**2
          CD(IJ) = ( XKAPPA / LOG( 1.0_JWRB + XNLEV/Z0ATM) )**2
        ENDDO
      ENDIF

IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)
! ----------------------------------------------------------------------
END SUBROUTINE OUTBETA
