SUBROUTINE HALPHAP(IJS, IJL, USTAR, UDIR, FL1, HALP)

! ----------------------------------------------------------------------

!**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER in the wind direction only


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

      USE YOWFRED  , ONLY : TH       , FR       , FR5      ,DELTH  ,FRIC     ,OLDWSFC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        , GM1      ,ZPI       ,ZPI4GM2,EPSMIN
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOWSHAL  , ONLY : TFAK     ,INDEP    ,CINV
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "femean.intfb.h"
#include "meansqs_lf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: HALP

      INTEGER(KIND=JWIM) :: IJ, K, M

      ! log of the conversion factor from windsea mean frequency to windsea peak frequency
      REAL(KIND=JWRB), PARAMETER :: XLOGMTP = LOG(0.85_JWRB)

      REAL(KIND=JWRB) :: CONST, COSPOS
      REAL(KIND=JWRB) :: COEF, WS, CHECKTA, XLOG
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ESEA, FSEA, XMSSSEA, FMAX, FPSEA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FLWS
!!debile
      REAL(KIND=JWRB) :: wage 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP_direct

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

!     COMPUTE THE PHILLIPS PARAMETER


!! Direct method:
      ALPHAP(:) = 0.0_JWRB
      CONST = DELTH*ZPI4GM2*FR5(NFRE)
      DO K = 1, NANG
        DO IJ = IJS, IJL
!!!          ALPHAP(IJ) = ALPHAP(IJ) + CONST*FL1(IJ,K,NFRE)
!         only in the wind sector
          COSPOS = 0.5_JWRB + SIGN(0.5_JWRB, COS(TH(K)-UDIR(IJ)) )
          ALPHAP(IJ) = ALPHAP(IJ) + COSPOS*CONST*FL1(IJ,K,NFRE)
        ENDDO
      ENDDO
!!debile
    alphap_direct(:) = alphap(:)
      

!! Via the mean square slope and mean frequency of the windsea:

      COEF = OLDWSFC*FRIC

      DO M = 1, NFRE
        DO IJ = IJS, IJL
          XINVWVAGE(IJ,M)=USTAR(IJ)*CINV(INDEP(IJ),M)
        ENDDO
      ENDDO

      DO K = 1, NANG
        DO IJ = IJS, IJL
          DIRCOEF(IJ,K)=COEF*COS(TH(K)-UDIR(IJ))
        ENDDO
      ENDDO

!     Find windsea spectrum
      FMAX(:) = 0.0_JWRB
      FPSEA(:) = FR(NFRE-1)
      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            WS = 0.5_JWRB + SIGN(0.5_JWRB, (CHECKTA-1.0_JWRB) )
            FLWS(IJ,K,M) = WS*FL1(IJ,K,M)
            IF(FLWS(IJ,K,M) > FMAX(IJ)) THEN
              FMAX(IJ) = FLWS(IJ,K,M)
              FPSEA(IJ) = FR(M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!!! totally crazy
!! Direct method:
      ALPHAP_direct(:) = 0.0_JWRB
      CONST = DELTH*ZPI4GM2*FR5(NFRE)
      DO K = 1, NANG
        DO IJ = IJS, IJL
          ALPHAP_direct(IJ) = ALPHAP_direct(IJ) + CONST*FLWS(IJ,K,NFRE)
        ENDDO
      ENDDO
!!debile
      


      CALL FEMEAN(FLWS, IJS, IJL, ESEA, FSEA)

      CALL MEANSQS_LF (NFRE, IJS, IJL, FLWS, XMSSSEA)

      ALPHAP(:) = 0.0_JWRB
      XLOG = LOG(FR(NFRE)) - XLOGMTP
      DO IJ = IJS, IJL
         IF(ESEA(IJ) > 0.0_JWRB ) THEN
           ALPHAP(IJ) = XMSSSEA(IJ) / ( XLOG - LOG(MIN(FSEA(IJ),FR(NFRE))) )
         ELSE
           ALPHAP(IJ) = 0.0065_JWRB
         ENDIF
!!debile
      wage = (g/(ZPI*0.85_JWRB*FSEA(IJ)))/ustar(IJ)
      write(*,*) 'debile halphap ', wage, alphap_direct(ij), ALPHAP(IJ), FSEA(IJ),FPSEA(IJ),ustar(ij)

!!!!!!!
      ENDDO


!    1/2 ALPHAP:
      HALP(:) = 0.5_JWRB*MIN(ALPHAP(:), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
