      SUBROUTINE STHQ (GFL, IJS, IJL, KIJS, KIJL, THQ)

! ----------------------------------------------------------------------

!**** *STHQ* - COMPUTATION OF MEAN WAVE DIRECTION AT EACH GRID POINT.

!     S.D. HASSELMANN
!     OPTIMIZED BY L. ZAMBRESKY

!*    PURPOSE.
!     --------

!       TO COMPUTE MEAN WAVE DIRECTION AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *STHQ (GFL, IJS, IJL, KIJS, KIJL, THQ)*
!          *GFL*    - SPECTRUM.
!          *IJS:IJL*- 1st DIMENSION of GFL
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *THQ*    - MEAN WAVE DIRECTION

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

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: THQ

      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, SI, CI

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STHQ',0,ZHOOK_HANDLE)

!*    1. INITIALISE SIN AND COS ARRAYS.
!        ------------------------------

      DO IJ=KIJS,KIJL
        SI(IJ) = 0.0_JWRB
        CI(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+GFL(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          SI(IJ) = SI(IJ)+SINTH(K)*TEMP(IJ)
          CI(IJ) = CI(IJ)+COSTH(K)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTE MEAN DIRECTION.
!        -----------------------

      DO IJ=KIJS,KIJL
        IF (CI(IJ).EQ.0.0_JWRB) CI(IJ) = EPSMIN
      ENDDO
      DO IJ=KIJS,KIJL
        THQ(IJ) = ATAN2(SI(IJ),CI(IJ))
      ENDDO
      DO IJ=KIJS,KIJL
        IF (THQ(IJ).LT.0.0_JWRB) THQ(IJ) = THQ(IJ) + ZPI
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STHQ',1,ZHOOK_HANDLE)

      END SUBROUTINE STHQ
