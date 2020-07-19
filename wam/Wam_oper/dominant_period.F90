      SUBROUTINE DOMINANT_PERIOD (F, IJS, IJL, DP)

! ----------------------------------------------------------------------

!**** *DOMINANT_PERIOD* - COMPUTATION OF THE DOMINANT PERIOD

!*    PURPOSE.
!     --------

!       COMPUTE DOMINANT PERIOD AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *DOMINANT_PERIOD (F, IJS, IJL, DP)*
!              *F*   - 2D-SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *DP*  - DOMINANT PERIOD

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

      USE YOWFRED  , ONLY : FR       ,DFIM_SIM ,DFIMFR_SIM,DELTH  ,WETAIL, WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: DP

      REAL(KIND=JWRB), PARAMETER :: FLTHRS = 0.1_JWRB

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP, EM, FMAX

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        EM(IJ) = 0.0_JWRB
        DP(IJ) = 0.0_JWRB
        FMAX(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF( F(IJ,K,M) > FMAX(IJ) ) THEN
               FMAX(IJ) = F(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      FMAX(:) = FLTHRS*FMAX(:)

      DO M=1,NFRE_ODD
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF( F(IJ,K,M) > FMAX(IJ) ) THEN
              TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          TEMP(IJ) = TEMP(IJ)**4
          EM(IJ) = EM(IJ)+DFIM_SIM(M)*TEMP(IJ)
          DP(IJ) = DP(IJ)+DFIMFR_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF(EM(IJ).GT.0.0_JWRB .AND. DP(IJ).GT.EPSMIN ) THEN
          DP(IJ) = EM(IJ)/DP(IJ)
        ELSE
          DP(IJ) = 0.0_JWRB 
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',1,ZHOOK_HANDLE)

      END SUBROUTINE DOMINANT_PERIOD
