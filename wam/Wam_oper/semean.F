      SUBROUTINE SEMEAN (F3, IJS, IJL, EM, LLEPSMIN)

! ----------------------------------------------------------------------

!**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.

!     S.D. HASSELMANN.
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEMEAN(F3, IJS, IJL, EM)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *EM*  - MEAN ENERGY
!          *LLEPSMIN* - TRUE IF THE WAVE ENERGY IS AT LEAST = EPSMIN

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,WETAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER :: IJ, M, K

      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: EM
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3

      REAL :: DELT25
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: TEMP

      LOGICAL, INTENT(IN) :: LLEPSMIN

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEMEAN',0,ZHOOK_HANDLE)
#endif


!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      IF(LLEPSMIN) THEN
        DO IJ=IJS,IJL
          EM(IJ) = EPSMIN
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          EM(IJ) = 0.
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,NFRE
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
          EM(IJ) = EM(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. ADD TAIL ENERGY.
!        ----------------

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DO IJ=IJS,IJL
        EM(IJ) = EM(IJ)+DELT25*TEMP(IJ)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEMEAN',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SEMEAN
