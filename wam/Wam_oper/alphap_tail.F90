SUBROUTINE ALPHAP_TAIL(IJS, IJL, FL1, ALPHAP)

! ----------------------------------------------------------------------

!**** *ALPHAP_TAIL* - COMPUTATION OF PHILLIPS PARAMETER using the last frequency bins


!**   INTERFACE.
!     ----------

!       *CALL* *ALPHAP_TAIL(IJS, IJL, FL1, ALPHAP)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *FL1*  - SPECTRA
!          *ALPHAP*  - PHILLIPS PARAMETER 

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH       ,FR5      ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: ALPHAP

      INTEGER(KIND=JWIM) :: IJ, K

      REAL(KIND=JWRB) :: CONST
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ALPHAP_TAIL',0,ZHOOK_HANDLE)

!     COMPUTE THE PHILLIPS PARAMETER
      CONST = DELTH*ZPI4GM2*FR5(NFRE)
      ALPHAP(:) = 0.0_JWRB
      DO K = 1, NANG
        DO IJ = IJS, IJL
          ALPHAP(IJ) = ALPHAP(IJ) + CONST*FL1(IJ,K,NFRE)
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('ALPHAP_TAIL',1,ZHOOK_HANDLE)

END SUBROUTINE ALPHAP_TAIL
