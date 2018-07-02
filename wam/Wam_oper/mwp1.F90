      SUBROUTINE MWP1 (F, IJS, IJL, EMEAN, MEANWP1)

! ----------------------------------------------------------------------

!**** *MWP1* - COMPUTATION OF MEAN PERIOD BASED ON THE FIRST MOMENT
!              FOR EACH GRID POINT.

!     J-R BIDLOT    ECMWF     MARCH 2000

!*    PURPOSE.
!     --------

!       COMPUTE MEAN PERIOD 1 AT EACH GRID POINT.
!       THE MEAN ENERGY MUST BE COMPUTED BEFORE CALLING THIS ROUTINE.

!**   INTERFACE.
!     ----------

!       *CALL* *MWP1 (F, IJS, IJL, EMEAN, MEANWP1)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EMEAN* - MEAN WAVE ENERGY FOR THE WAVE SYSTEM (INPUT).
!              *MEANWP1* - MEAN PERIOD BASED ON THE FIRST MOMENT.

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

      USE YOWFRED  , ONLY : FR       ,DFIMFR_SIM,DELTH    ,WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: MEANWP1

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: COEF_FR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MWP1',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        MEANWP1(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          MEANWP1(IJ) = MEANWP1(IJ)+DFIMFR_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      COEF_FR = WP1TAIL*DELTH*FR(NFRE_ODD)**2

      DO IJ=IJS,IJL
        MEANWP1(IJ) = MEANWP1(IJ)+COEF_FR*TEMP(IJ)
        IF(EMEAN(IJ).GT.EPSMIN .AND. MEANWP1(IJ).GT.EPSMIN ) THEN
          MEANWP1(IJ) = EMEAN(IJ)/MEANWP1(IJ)
        ELSE
          MEANWP1(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MWP1',1,ZHOOK_HANDLE)

      END SUBROUTINE MWP1
