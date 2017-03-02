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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFFR  ,DELTH    ,WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,K,M,IJS,IJL
      REAL :: F(IJS:IJL,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL) :: TEMP
      REAL,DIMENSION(IJS:IJL), intent(in) :: EMEAN
      REAL,DIMENSION(IJS:IJL), intent(out) :: MEANWP1
      REAL :: COEF_FR

! ----------------------------------------------------------------------

      DO IJ=IJS,IJL
        MEANWP1(IJ) = 0.
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          MEANWP1(IJ) = MEANWP1(IJ)+DFFR(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      COEF_FR = WP1TAIL*DELTH*FR(NFRE)**2

      DO IJ=IJS,IJL
        MEANWP1(IJ) = MEANWP1(IJ)+COEF_FR*TEMP(IJ)
        IF(EMEAN(IJ).GT.EPSMIN .AND. MEANWP1(IJ).GT.EPSMIN ) THEN
          MEANWP1(IJ) = EMEAN(IJ)/MEANWP1(IJ)
        ELSE
          MEANWP1(IJ) = 0. 
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE MWP1
