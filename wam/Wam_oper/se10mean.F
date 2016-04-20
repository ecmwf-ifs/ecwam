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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER :: IJ, M, K, MCUT

      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: E10
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3

      REAL, PARAMETER :: FCUT=0.1
      REAL :: DFCUT
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NFRE) :: DFIMLOC
      REAL, DIMENSION(IJS:IJL) :: TEMP


! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SE10MEAN',0,ZHOOK_HANDLE)
#endif

      MCUT=NINT(LOG10(FCUT/FR(1))*FLOGSPRDM1)+1
      MCUT=MIN(MAX(1,MCUT),NFRE)

      DFIMLOC(:)=DFIM(:)
      DFCUT=0.5*MAX(0.,FCUT-0.5*(FR(MAX(1,MCUT-1))+FR(MCUT)))
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

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SE10MEAN',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SE10MEAN
