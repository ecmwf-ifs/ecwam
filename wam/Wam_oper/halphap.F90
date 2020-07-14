SUBROUTINE HALPHAP(IJS, IJL, FL1, THW, HALP)

! ----------------------------------------------------------------------

!**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER


!**   INTERFACE.
!     ----------

!       *CALL* *HALPHAP(IJS, IJL, FL1, THW, HALP)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *FL1*  - SPECTRA
!          *THW*  - WIND DIRECTION
!          *HALP*   - 1/2 PHILLIPS PARAMETER 

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH       , FR5      ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        , ZPI      ,ZPI4GM2
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: HALP

      REAL(KIND=JWRB), PARAMETER :: ALPHAPMIN=0.004_JWRB
      INTEGER(KIND=JWIM) :: IJ, K
      INTEGER(KIND=JWIM) :: IFRPH

      REAL(KIND=JWRB) :: CONST
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

!     COMPUTE THE PHILLIPS PARAMETER (only in the wind direction)
      IFRPH=NFRE
      CONST = DELTH*ZPI4GM2*FR5(IFRPH)
      ALPHAP(:) = 0.0_JWRB
      DO K = 1, NANG
        DO IJ = IJS, IJL
          ALPHAP(IJ) = ALPHAP(IJ) + CONST*FL1(IJ,K,IFRPH)*SIGN(1.0_JWRB, COS(TH(K)-THW(IJ)) )
        ENDDO
      ENDDO

!!!! debile
      write(*,*) 'alphap ', ALPHAP(1)



      HALP(:) = 0.5_JWRB*MIN(MAX(ALPHAP(:),ALPHAPMIN), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
