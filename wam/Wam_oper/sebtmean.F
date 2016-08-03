      SUBROUTINE SEBTMEAN (F3, IJS, IJL, TB, TT, EBT)

! ----------------------------------------------------------------------

!**** *SEBTMEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                  WITH  TB <= PERIODS <= TT
!                  NO HIGH FREQUENCY TAIL IS ASSUMED !!!!

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEBTMEAN(F3, IJS, IJL, TB, TT, EBT)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *TB*  - BOTTOM PERIOD
!          *TT*  - TOP PERIOD
!          *EBT* - MEAN ENERGY

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, TB, TT
      INTEGER :: IJ, M, K, MCUTB, MCUTT

      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: EBT
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3

      REAL :: FCUTB, FCUTT
      REAL :: DFCUT
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NFRE) :: DFIMLOC
      REAL, DIMENSION(IJS:IJL) :: TEMP


! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',0,ZHOOK_HANDLE)
#endif

      DFIMLOC(:)=DFIM(:)

      FCUTB=MAX(FR(1),MIN(1./MAX(TB,EPSMIN),FR(NFRE)))
      FCUTT=MAX(FR(1),MIN(1./MAX(TT,EPSMIN),FR(NFRE)))

      MCUTB=NINT(LOG10(FCUTB/FR(1))*FLOGSPRDM1)+1
      MCUTB=MIN(MAX(1,MCUTB),NFRE)
      DFCUT=0.5*MAX(0.,FCUTB-0.5*(FR(MAX(1,MCUTB-1))+FR(MCUTB)))
      DFIMLOC(MCUTB)=MIN(DFCUT*DELTH,DFIM(MCUTB))
      DFIMLOC(MCUTB)=DFIM(MCUTB)-DFIMLOC(MCUTB)

      MCUTT=NINT(LOG10(FCUTT/FR(1))*FLOGSPRDM1)+1
      MCUTT=MIN(MAX(1,MCUTT),NFRE)
      DFCUT=0.5*MAX(0.,FCUTT-0.5*(FR(MAX(1,MCUTT-1))+FR(MCUTT)))
      DFIMLOC(MCUTT)=MIN(DFCUT*DELTH,DFIM(MCUTT))

      IF(FCUTB.EQ.FCUTT) MCUTT=MCUTB-1

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=IJS,IJL
        EBT(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=MCUTB,MCUTT
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
          EBT(IJ) = EBT(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SEBTMEAN
