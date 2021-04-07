      SUBROUTINE MWP2 (F, IJS, IJL, MEANWP2)

! ----------------------------------------------------------------------

!**** *MWP2* - COMPUTATION OF MEAN PERIOD BASED ON THE SECOND MOMENT
!              FOR EACH GRID POINT.

!     J-R BIDLOT    ECMWF     MARCH 2000

!*    PURPOSE.
!     --------

!       COMPUTE MEAN PERIOD 2 AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MWP2 (F, IJS, IJL, MEANWP2)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *MEANWP2* - MEAN PERIOD BASED ON THE SECOND MOMENT.

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

      USE YOWFRED  , ONLY : FR       ,DFIM_SIM ,DFIMFR2_SIM,DELTH  ,WETAIL, WP2TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: MEANWP2

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELT25, COEF_FR, FR1M1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP, EM

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MWP2',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        EM(IJ) = 0.0_JWRB
        MEANWP2(IJ) = 0.0_JWRB
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
          EM(IJ) = EM(IJ)+DFIM_SIM(M)*TEMP(IJ)
          MEANWP2(IJ) = MEANWP2(IJ)+DFIMFR2_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      FR1M1 = 1.0_JWRB/FR(1)
      DELT25 = WETAIL*FR(NFRE_ODD)*DELTH
      COEF_FR = WP2TAIL*DELTH*FR(NFRE_ODD)**3

      DO IJ=IJS,IJL
        EM(IJ) = EM(IJ)+DELT25*TEMP(IJ)
        MEANWP2(IJ) = MEANWP2(IJ)+COEF_FR*TEMP(IJ)
        IF(EM(IJ).GT.0.0_JWRB .AND. MEANWP2(IJ).GT.EPSMIN ) THEN
          MEANWP2(IJ) = SQRT(EM(IJ)/MEANWP2(IJ))
          MEANWP2(IJ) = MIN(MEANWP2(IJ),FR1M1)
        ELSE
          MEANWP2(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MWP2',1,ZHOOK_HANDLE)

      END SUBROUTINE MWP2
