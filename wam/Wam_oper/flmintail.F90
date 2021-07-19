SUBROUTINE FLMINTAIL(IJS, IJL, U10, UDIR, USTAR, FMEANWS, FL1)

! ----------------------------------------------------------------------

!**** *FLMINTAIL* - SETS MINIMUM VALUE FOR THE LAST FREQUENCY SPECTRAL BIN.


!**   INTERFACE.
!     ----------

!       *CALL* *FLMINTAIL(IJS, IJL, UDIR, USTAR, FMEANWS,FL1)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *U10* - WIND SPEED
!          *UDIR* - WIND SPEED DIRECTION
!          *USTAR*- FRICTION VELOCITY
!          *FMEANWS* - MEAN FREQUENCY OF WINDSEA
!          *FL1*  - SPECTRA

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPHYS  , ONLY : FLMINFAC 
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: FMEANWS 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM) :: IJ, K

      REAL(KIND=JWRB) :: COSPOS, COS2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FMIN

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FLMINTAIL',0,ZHOOK_HANDLE)

      DO IJ = IJS, IJL
        IF(U10(IJ) > 5.0_JWRB) THEN
          FMIN(IJ) = FLMINFAC * FMEANWS(IJ) * USTAR(IJ)
        ELSE
!          FMIN(IJ) = 0.0_JWRB 
          FMIN(IJ) = FLMINFAC * FMEANWS(IJ) * USTAR(IJ)
        ENDIF
      ENDDO

      DO K = 1, NANG
        DO IJ = IJS, IJL
!!          COSPOS = 0.5_JWRB + SIGN(0.5_JWRB, COS(TH(K)-UDIR(IJ)) )
!!          FL1(IJ,K,NFRE) = MAX(FL1(IJ,K,NFRE),FMIN(IJ)*COSPOS) 
          COS2 = MAX(0.0_JWRB,COS(TH(K)-UDIR(IJ)))**2
          FL1(IJ,K,NFRE) = MAX(FL1(IJ,K,NFRE),FMIN(IJ)*COS2) 
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('FLMINTAIL',1,ZHOOK_HANDLE)

END SUBROUTINE FLMINTAIL
