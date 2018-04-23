      SUBROUTINE FEMEAN (F, IJS, IJL, EM, FM)

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

!       *CALL* *FEMEAN (F, IJS, IJL, EM, FM)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY (OUTPUT)
!              *FM*  - MEAN WAVE FREQUENCY (OUTPUT)

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
     &                WETAIL    ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: EM, FM

      INTEGER :: IJ, M, K
      REAL :: DELT25, DELT2, DEL2
      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(IJS:IJL) :: TEMP2

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FEMEAN',0,ZHOOK_HANDLE)
#endif

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = 0.
        FM(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        K=1
        DO IJ=IJS,IJL
          TEMP2(IJ) = F(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
          FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
        ENDDO
      ENDDO


!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
        FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
        FM(IJ) = EM(IJ)/FM(IJ)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FEMEAN',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE FEMEAN
