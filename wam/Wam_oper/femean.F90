      SUBROUTINE FEMEAN (F, IJS, IJL, KIJS, KIJL, EM, FM)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF VARIANCE and MEAN FREQUENCY AT EACH GRID POINT

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL, KIJS, KIJL, EM, FM)*
!              *F*       - SPECTRUM.
!              *IJS:IJL* - 1st DIMENSION of F
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *EM*      - MEAN WAVE ENERGY (OUTPUT)
!              *FM*      - MEAN WAVE FREQUENCY (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,    &
     &                WETAIL    ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EM, FM


      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JWRB) :: DELT25, DELT2, DEL2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FEMEAN',0,ZHOOK_HANDLE)

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=KIJS,KIJL
        EM(IJ) = 0.0_JWRB
        FM(IJ) = 0.0_JWRB
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        K=1
        DO IJ=KIJS,KIJL
          TEMP2(IJ) = MAX(F(IJ,K,M), EPSMIN)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP2(IJ) = TEMP2(IJ)+ MAX(F(IJ,K,M), EPSMIN)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
          FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
        ENDDO
      ENDDO


!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=KIJS,KIJL
        EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
        FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
        FM(IJ) = EM(IJ)/FM(IJ)
        FM(IJ) = MAX(FM(IJ),FR(1))
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FEMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE FEMEAN
