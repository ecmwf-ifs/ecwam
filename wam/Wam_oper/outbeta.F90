SUBROUTINE OUTBETA (KIJS, KIJL, PRCHAR, FF_NOW, BETAHQ)

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

!       *CALL* *OUTBETA (KIJS, KIJL, PRCHAR, FF_NOW, BETAHQ)
!         *KIJS*    - INDEX OF FIRST GRIDPOINT.
!         *KIJL*    - INDEX OF LAST GRIDPOINT.
!         *PRCHAR* - DEFAULT VALUE FOR CHARNOCK
!         *FF_NOW* - FORCING F
!         *U10*    - WIND SPEED IN M/S.
!         *US*     - FRICTION VELOCITY IN M/S.
!         *Z0*     - ROUGHNESS LENGTH IN M.
!         *Z0B*    - BACKGROUND ROUGHNESS LENGTH IN M.
!         *BETA*   - CHARNOCK FIELD (BLOCK ARRAY)
!         *BETAHQ* - EQUIVALENT CHARNOCK FIELD FOR HEAT AND MOISTURE
!                    (i.e. Charnock with the background roughness removed)   
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
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOWPHYS  , ONLY : RNUM     ,ALPHAMIN  , ALPHAMAX
      USE YOWICE   , ONLY : LICERUN  ,LWAMRSETCI, CITHRSH

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: PRCHAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: FF_NOW

      REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
      REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB
      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: Z0VIS, ZN, USM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL)  :: ALPHAMAXU10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

ASSOCIATE(U10 => FF_NOW%WSWAVE, &
 &        US => FF_NOW%UFRIC, &
 &        Z0 => FF_NOW%Z0M, &
 &        Z0B => FF_NOW%Z0B, &
 &        BETA => FF_NOW%CHNK )


!*    COMPUTE CHARNOCK 'CONSTANT' BETA.
!     ---------------------------------

      IF (LLGCBZ0) THEN
        ZN = RNUM
        ALPHAMAXU10(:)=ALPHAMAX
      ELSE
        ZN = 0.0_JWRB
        ALPHAMAXU10(:)=MIN(ALPHAMAX,AMAX+BMAX*U10(:))
      ENDIF

      DO IJ = KIJS,KIJL
        USM = 1.0_JWRB/MAX(US(IJ),EPSUS)
        Z0VIS = ZN*USM
        BETA(IJ) = G*(Z0(IJ)-Z0VIS)*USM**2
        BETA(IJ) = MAX(MIN(BETA(IJ),ALPHAMAXU10(IJ)),ALPHAMIN)
        BETAHQ(IJ) = MAX(BETA(IJ)-G*Z0B(IJ)*USM**2,ALPHAMIN)
      ENDDO

END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)
! ----------------------------------------------------------------------
END SUBROUTINE OUTBETA
