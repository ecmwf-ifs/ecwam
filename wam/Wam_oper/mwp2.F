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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFFR2   ,DELTH  ,WP2TAIL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,K,M,IJS,IJL
      REAL :: F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL) :: TEMP, MEANWP2, EMEAN
      REAL :: COEF_FR

! ----------------------------------------------------------------------

      DO IJ=IJS,IJL
        MEANWP2(IJ) = 0.
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
          MEANWP2(IJ) = MEANWP2(IJ)+DFFR2(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      COEF_FR = WP2TAIL*DELTH*FR(NFRE)**3

      DO IJ=IJS,IJL
        MEANWP2(IJ) = MEANWP2(IJ)+COEF_FR*TEMP(IJ)
        IF(EMEAN(IJ).GT.EPSMIN .AND. MEANWP2(IJ).GT.EPSMIN ) THEN
          MEANWP2(IJ) = SQRT(EMEAN(IJ)/MEANWP2(IJ))
        ELSE
          MEANWP2(IJ) = 0. 
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE MWP2
