      SUBROUTINE DOMINANT_PERIOD (GFL, IJL, IJL, KIJL, KIJL, DP)

! ----------------------------------------------------------------------

!**** *DOMINANT_PERIOD* - COMPUTATION OF THE DOMINANT PERIOD

!*    PURPOSE.
!     --------

!       COMPUTE DOMINANT PERIOD AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *DOMINANT_PERIOD (GFL, IJL, IJL, KIJL, KIJL, DP)*
!              *GFL*     - 2D-SPECTRUM.
!              *IJS:IJL* - 1st DIMEMSION of GFL
!              *KIJL*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *DP*      - DOMINANT PERIOD

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMFR   ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJL, KIJL

      REAL(KIND=JWRB), DIMENSION(IJL:IJL,NANG,NFRE), INTENT(IN) :: GFL 

      REAL(KIND=JWRB), DIMENSION(KIJL:KIJL), INTENT(OUT) :: DP


      REAL(KIND=JWRB), PARAMETER :: FLTHRS = 0.1_JWRB

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL:KIJL) :: TEMP, EM, FCROP
      REAL(KIND=JWRB), DIMENSION(KIJL:KIJL,NFRE) :: F1D4

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',0,ZHOOK_HANDLE)

      F1D4(:,:) = 0.0_JWRB
      EM(:) = 0.0_JWRB
      DP(:) = 0.0_JWRB
      FCROP(:) = 0.0_JWRB

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJL,KIJL
            IF( GFL(IJ,K,M) > FCROP(IJ) ) THEN
               FCROP(IJ) = GFL(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      FCROP(:) = FLTHRS*FCROP(:)

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJL,KIJL
            IF(GFL(IJ,K,M) > FCROP(IJ)) F1D4(IJ,M) = F1D4(IJ,M)+GFL(IJ,K,M)*DELTH
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJL,KIJL
          F1D4(IJ,M) = F1D4(IJ,M)**4
          EM(IJ) = EM(IJ)+DFIM(M)*F1D4(IJ,M)
          DP(IJ) = DP(IJ)+DFIMFR(M)*F1D4(IJ,M)
        ENDDO
      ENDDO

      DO IJ=KIJL,KIJL
        IF(EM(IJ) > 0.0_JWRB .AND. DP(IJ).GT.EPSMIN ) THEN
          DP(IJ) = EM(IJ)/DP(IJ)
        ELSE
          DP(IJ) = 0.0_JWRB 
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('DOMINANT_PERIOD',1,ZHOOK_HANDLE)

      END SUBROUTINE DOMINANT_PERIOD
