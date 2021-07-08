      SUBROUTINE PEAKFRI (F, IJS, IJL, IPEAKF, EPEAKF, F1D)

! ----------------------------------------------------------------------

!**** *PEAKFRI* - COMPUTATION OF INDEX OF THE PEAK FREQUENCY AND
!                 THE CORRESPONDING VALUE OF THE 1-D SPECTRAM
!                 AT EACH GRID POINT.

!     JEAN BIDLOT    ECMWF           MARCH 2000 

!*    PURPOSE.
!     --------

!       COMPUTE PEAK FREQUENCY INDEX AT EACH GRID POINT.
!       IF NOT PEAK FOUND THAN EPEAKF WILL BE SET TO 0. AND
!       IPEAKF TO NFRE

!**   INTERFACE.
!     ----------

!       *CALL* *PEAKFRI (F, IJS, IJL, IPEAKF, EPEAKF, F1D)*
!          *F*     - SPECTRUM.
!          *IJS*   - INDEX OF FIRST GRIDPOINT
!          *IJL*   - INDEX OF LAST GRIDPOINT
!          *IPEAKF - INDEX OF PEAK FREQUENCY 
!          *EPEAKF* - 1-D SPECTRUM AT PEAK FREQUENCY
!          *F1D*   - 1D FREQUENCY SPECTRUM

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

      USE YOWFRED  , ONLY : DELTH    ,FR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: IPEAKF(IJS:IJL) 
      REAL(KIND=JWRB), INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: EPEAKF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: F1D

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PEAKFRI',0,ZHOOK_HANDLE)

!*    1. INITIALIZE ARRAYS
!        -----------------

      DO IJ = IJS,IJL
        EPEAKF(IJ) = 0.0_JWRB
        IPEAKF(IJ) = NFRE 
      ENDDO

!*    2. LOOP OVER FREQUENCIES
!        ---------------------

      DO M = 1, NFRE

!*    2.1 COMPUTE 1-D SPECTRUM
!        ---------------------

        DO IJ = IJS,IJL
          F1D(IJ,M) = 0.0_JWRB
        ENDDO
        DO K = 1,NANG
          DO IJ = IJS,IJL
            F1D(IJ,M) = F1D(IJ,M) + F(IJ,K,M)*DELTH
          ENDDO
        ENDDO

!*    2.2 DEFINE PEAK FREQUENCY
!         ---------------------

        DO IJ = IJS,IJL
          IF (EPEAKF(IJ).LT.F1D(IJ,M)) THEN
            EPEAKF(IJ) = F1D(IJ,M)
            IPEAKF(IJ) = M 
          ENDIF
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('PEAKFRI',1,ZHOOK_HANDLE)

      END SUBROUTINE PEAKFRI
