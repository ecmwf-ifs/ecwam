      SUBROUTINE MWP2 (F, IJS, IJL, EMEAN, MEANWP2)

! ----------------------------------------------------------------------

!**** *MWP2* - COMPUTATION OF MEAN PERIOD BASED ON THE SECOND MOMENT
!              FOR EACH GRID POINT.

!     J-R BIDLOT    ECMWF     MARCH 2000

!*    PURPOSE.
!     --------

!       COMPUTE MEAN PERIOD 2 AT EACH GRID POINT.
!       THE MEAN ENERGY MUST BE COMPUTED BEFORE CALLING THIS ROUTINE.

!**   INTERFACE.
!     ----------

!       *CALL* *MWP2 (F, IJS, IJL, EMEAN, MEANWP2)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EMEAN* - MEAN WAVE ENERGY FOR THE WAVE SYSTEM (INPUT).
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

      USE YOWFRED  , ONLY : FR       ,DFIMFR2_SIM,DELTH  ,WP2TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: MEANWP2

      INTEGER :: IJ, K, M
      REAL,DIMENSION(IJS:IJL) :: TEMP
      REAL :: ZHOOK_HANDLE
      REAL :: COEF_FR

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MWP2',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        MEANWP2(IJ) = 0.
      ENDDO

      DO M=1,NFRE_ODD
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          MEANWP2(IJ) = MEANWP2(IJ)+DFIMFR2_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      COEF_FR = WP2TAIL*DELTH*FR(NFRE_ODD)**3

      DO IJ=IJS,IJL
        MEANWP2(IJ) = MEANWP2(IJ)+COEF_FR*TEMP(IJ)
        IF(EMEAN(IJ).GT.EPSMIN .AND. MEANWP2(IJ).GT.EPSMIN ) THEN
          MEANWP2(IJ) = SQRT(EMEAN(IJ)/MEANWP2(IJ))
        ELSE
          MEANWP2(IJ) = 0. 
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MWP2',1,ZHOOK_HANDLE)
      END SUBROUTINE MWP2
