      SUBROUTINE SE10MEAN (F3, IJS, IJL, E10)

! ----------------------------------------------------------------------

!**** *SE10MEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                 WITH PERIOD LARGER THAN 10s

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SE10MEAN(F3, IJS, IJL, E10)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *E10* - MEAN ENERGY

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: F3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: E10


      INTEGER(KIND=JWIM) :: IJ, M, K, MCUT

      REAL(KIND=JWRB), PARAMETER :: FCUT=0.1_JWRB
      REAL(KIND=JWRB) :: DFCUT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFIMLOC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',0,ZHOOK_HANDLE)

      MCUT=NINT(LOG10(FCUT/FR(1))*FLOGSPRDM1)+1
      MCUT=MIN(MAX(1,MCUT),NFRE)

      DFIMLOC(:)=DFIM(:)
      DFCUT=0.5_JWRB*MAX(0.0_JWRB,FCUT-0.5_JWRB*(FR(MAX(1,MCUT-1))+FR(MCUT)))
      DFIMLOC(MCUT)=MIN(DFCUT*DELTH,DFIM(MCUT))


!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=IJS,IJL
        E10(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,MCUT
        K=1
        DO IJ=IJS,IJL
          TEMP(IJ) = F3(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          E10(IJ) = E10(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SE10MEAN
