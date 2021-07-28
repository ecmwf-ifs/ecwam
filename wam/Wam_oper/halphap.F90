SUBROUTINE HALPHAP(IJS, IJL, USTAR, UDIR, FL1, HALP)

! ----------------------------------------------------------------------

!**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER


!**   INTERFACE.
!     ----------

!       *CALL* *HALPHAP(IJS, IJL, UDIR, FL1, HALP)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *UDIR* - WIND SPEED DIRECTION
!          *USTAR* - FRICTION VELOCITY
!          *FL1*  - SPECTRA
!          *HALP*   - 1/2 PHILLIPS PARAMETER 

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       , TH       , FR5      ,DELTH
      USE YOWPARAM , ONLY : NANG     , NFRE
      USE YOWPCONS , ONLY : G        , ZPI      ,ZPI4GM2
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "femean.intfb.h"
#include "meansqs_lf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
!!! debile: it it works remove use of ustar and udir
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: HALP

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: F1D
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XMSS, EM, FM 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: WD 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FLWD

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

      ! Find spectrum in wind direction
      DO K = 1, NANG
        DO IJ = IJS, IJL
           WD(IJ,K) = 0.5_JWRB + 0.5_JWRB * SIGN(1.0_JWRB, COS(TH(K)-UDIR(IJ)) )
        ENDDO
      ENDDO

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
             FLWD(IJ,K,M) = FL1(IJ,K,M) * WD(IJ,K) 
          ENDDO
        ENDDO
      ENDDO

      CALL MEANSQS_LF(NFRE, IJS, IJL, FLWD, XMSS)
      CALL FEMEAN (FLWD, IJS, IJL, EM, FM)

      DO IJ = IJS, IJL
        IF(EM(IJ) > 0.0_JWRB .AND. FM(IJ) < FR(NFRE-2) ) THEN
          ALPHAP(IJ) = XMSS(IJ) /LOG(FR(NFRE)/FM(IJ))
          IF ( ALPHAP(IJ) > ALPHAPMAX ) THEN
            ! some odd cases, revert to tail value
            F1D = 0.0_JWRB
            DO K = 1, NANG
              F1D = F1D + FLWD(IJ,K,NFRE)*DELTH
            ENDDO
            ALPHAP(IJ) = ZPI4GM2*FR5(NFRE)*F1D
          ENDIF
        ELSE
          F1D = 0.0_JWRB
          DO K = 1, NANG
            F1D = F1D + FLWD(IJ,K,NFRE)*DELTH
          ENDDO
          ALPHAP(IJ) = ZPI4GM2*FR5(NFRE)*F1D
        ENDIF
      ENDDO

!     1/2 ALPHAP:
      HALP(:) = 0.5_JWRB*MIN(ALPHAP(:), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
