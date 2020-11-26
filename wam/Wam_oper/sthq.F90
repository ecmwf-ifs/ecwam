      SUBROUTINE STHQ (F3, IJS, IJL, THQ)

! ----------------------------------------------------------------------

!**** *STHQ* - COMPUTATION OF MEAN WAVE DIRECTION AT EACH GRID POINT.

!     S.D. HASSELMANN
!     OPTIMIZED BY L. ZAMBRESKY

!*    PURPOSE.
!     --------

!       TO COMPUTE MEAN WAVE DIRECTION AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *STHQ (F3, IJS, IJL, THQ)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *THQ* - MEAN WAVE DIRECTION

!     METHOD.
!     -------

!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,EPSMIN
      USE YOWFRED  , ONLY : DFIM     ,COSTH    ,SINTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: THQ

      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP, SI, CI

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STHQ',0,ZHOOK_HANDLE)

!*    1. INITIALISE SIN AND COS ARRAYS.
!        ------------------------------

      DO IJ=IJS,IJL
        SI(IJ) = 0.0_JWRB
        CI(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          SI(IJ) = SI(IJ)+SINTH(K)*TEMP(IJ)
          CI(IJ) = CI(IJ)+COSTH(K)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTE MEAN DIRECTION.
!        -----------------------

      DO IJ=IJS,IJL
        IF (CI(IJ).EQ.0.0_JWRB) CI(IJ) = EPSMIN
      ENDDO
      DO IJ=IJS,IJL
        THQ(IJ) = ATAN2(SI(IJ),CI(IJ))
      ENDDO
      DO IJ=IJS,IJL
        IF (THQ(IJ).LT.0.0_JWRB) THQ(IJ) = THQ(IJ) + ZPI
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STHQ',1,ZHOOK_HANDLE)

      END SUBROUTINE STHQ
