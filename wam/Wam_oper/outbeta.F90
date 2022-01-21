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
!         *FF_NOW* - FORCING FIELDS
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
      USE YOWPHYS  , ONLY : RNUM     ,ALPHAMIN  , ALPHAMAX, ALPHA
      USE YOWICE   , ONLY : LICERUN  ,LWAMRSETCI, CITHRSH

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: PRCHAR
      TYPE(FORCING_FIELDS), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: FF_NOW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: BETAHQ 


      REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
      REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB

!     BETAHQ_REDUCE was introduced to reduce the impact of the sea state dependent
!     heat and moisture surface flux. It is just a dirty tuning knob !!
      REAL(KIND=JWRB), PARAMETER :: BETAHQ_REDUCE=0.9_JWRB

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: Z0VIS, ZN, USM, GUSM2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL)  :: ALPHAMAXU10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

ASSOCIATE(WSWAVE => FF_NOW%WSWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        Z0M => FF_NOW%Z0M, &
 &        Z0B => FF_NOW%Z0B, &
 &        CHNK => FF_NOW%CHNK )


!*    COMPUTE CHARNOCK 'CONSTANT' CHNK.
!     ---------------------------------

      IF (LLGCBZ0) THEN
        ZN = RNUM
        ALPHAMAXU10(:)=ALPHAMAX
      ELSE
        ZN = 0.0_JWRB
        ALPHAMAXU10(:)=MIN(ALPHAMAX,AMAX+BMAX*WSWAVE(:))
      ENDIF

      DO IJ = KIJS,KIJL
        USM = 1.0_JWRB/MAX(UFRIC(IJ), EPSUS)
        Z0VIS = ZN*USM
        GUSM2 = G*USM**2
        CHNK(IJ) = (Z0M(IJ)-Z0VIS)*GUSM2
        CHNK(IJ) = MAX(MIN(CHNK(IJ),ALPHAMAXU10(IJ)),ALPHAMIN)
        BETAHQ(IJ) = BETAHQ_REDUCE*MAX(CHNK(IJ)-MAX(Z0B(IJ)*GUSM2, ALPHA), ALPHAMIN)
      ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)
! ----------------------------------------------------------------------
END SUBROUTINE OUTBETA
