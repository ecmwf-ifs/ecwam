SUBROUTINE FLMINTAIL(IJS, IJL, KIJS, KIJL, U10, UDIR, USTAR, FMEANWS, GFL)

! ----------------------------------------------------------------------

!**** *FLMINTAIL* - SETS MINIMUM VALUE FOR THE LAST FREQUENCY SPECTRAL BIN.


!**   INTERFACE.
!     ----------

!       *CALL* *FLMINTAIL(KIJS, KIJL, UDIR, USTAR, FMEANWS,GFL)
!          *KIJS* - INDEX OF FIRST GRIDPOINT
!          *KIJL* - INDEX OF LAST GRIDPOINT
!          *U10* - WIND SPEED
!          *UDIR* - WIND SPEED DIRECTION
!          *USTAR*- FRICTION VELOCITY
!          *FMEANWS* - MEAN FREQUENCY OF WINDSEA
!          *GFL*  - SPECTRA

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

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: U10
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: FMEANWS 

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: GFL

      INTEGER(KIND=JWIM) :: IJ, K

      REAL(KIND=JWRB) :: COSPOS
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FMIN

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FLMINTAIL',0,ZHOOK_HANDLE)

      DO IJ = KIJS, KIJL
        IF(U10(IJ) > 5.0_JWRB) THEN
          FMIN(IJ) = FLMINFAC * FMEANWS(IJ) * USTAR(IJ)
        ELSE
          FMIN(IJ) = 0.0_JWRB 
        ENDIF
      ENDDO

      DO K = 1, NANG
        DO IJ = KIJS, KIJL
          COSPOS = 0.5_JWRB + SIGN(0.5_JWRB, COS(TH(K)-UDIR(IJ)) )
          GFL(IJ,K,NFRE) = MAX(GFL(IJ,K,NFRE),FMIN(IJ)*COSPOS) 
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('FLMINTAIL',1,ZHOOK_HANDLE)

END SUBROUTINE FLMINTAIL
