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

      ! conversion factor from windsea mean frequency to windsea peak frequency
      REAL(KIND=JWRB), PARAMETER :: XMTP = 0.85_JWRB
      REAL(KIND=JWRB), PARAMETER :: XLOGMTP = LOG(XMTP)
      ! Maximum value for the mss calculation: XMTPCUT* windsea mean frequency
      REAL(KIND=JWRB), PARAMETER :: XMTPCUT = 2.0_JWRB*XMTP
      REAL(KIND=JWRB), PARAMETER :: XLOGMTPCUT = LOG(XMTPCUT)

      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MMAX

      REAL(KIND=JWRB) :: CONST, COSPOS
      REAL(KIND=JWRB) :: XLOGFS, FCUT, XMSS_TAIL
      REAL(KIND=JWRB) :: COEF, WS, CHECKTA, XLOG, XLOGCUT
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ESEA, FSEA, XMSSSEA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FMAX, FPSEA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP_NFRE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FLWS
!!debile
      REAL(KIND=JWRB) :: wage 

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

!     COMPUTE THE PHILLIPS PARAMETER

      XLOGFS = LOG(FR(NFRE))

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
      MMAX(:) = NFRE
      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            WS = 0.5_JWRB + SIGN(0.5_JWRB, (CHECKTA-1.0_JWRB) )
            FLWS(IJ,K,M) = WS*FL1(IJ,K,M)
            IF(FLWS(IJ,K,M) > FMAX(IJ)) THEN
              FMAX(IJ) = FLWS(IJ,K,M)
              FPSEA(IJ) = FR(M)
              MMAX(IJ) = M
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!!!debile
    MMAX(:) = MIN(NFRE,MMAX(:)+8) 

!! Direct method based on last discretised spectral value
      ALPHAP_NFRE(:) = 0.0_JWRB
      DO K = 1, NANG
        DO IJ = IJS, IJL
          CONST = DELTH*ZPI4GM2*FR5(MMAX(IJ))
          ALPHAP_NFRE(IJ) = ALPHAP_NFRE(IJ) + CONST*FLWS(IJ,K,MMAX(IJ))
        ENDDO
      ENDDO


      CALL FEMEAN(FLWS, IJS, IJL, ESEA, FSEA)

      CALL MEANSQS_LF (NFRE, IJS, IJL, FLWS, XMSSSEA)


      ALPHAP(:) = 0.0_JWRB
      XLOG = XLOGFS - XLOGMTP
      XLOGCUT = XLOGMTPCUT - XLOGFS
      DO IJ = IJS, IJL
         IF(ESEA(IJ) > 0.0_JWRB ) THEN
           XMSS_TAIL = ALPHAP_NFRE(IJ) * MAX((LOG(FSEA(IJ)) + XLOGCUT), 0.0_JWRB)
           ALPHAP(IJ) = XMSS_TAIL + XMSSSEA(IJ) / ( XLOG - LOG(MIN(FSEA(IJ),FR(NFRE))) )
         ELSE
           ALPHAP(IJ) = 0.0065_JWRB
         ENDIF
!!debile
wage = (g/(ZPI*XMTP*FSEA(IJ)))/ustar(IJ)
write(*,*) 'debile halphap ', wage, ALPHAP_NFRE(ij), ALPHAP(IJ), FSEA(IJ),FPSEA(IJ),ustar(ij),mmax(ij)

!!!!!!!
      ENDDO


!    1/2 ALPHAP:
      HALP(:) = 0.5_JWRB*MIN(ALPHAP(:), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
