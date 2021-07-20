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

      USE YOWFRED  , ONLY : TH       , FR       , FR5      ,FRATIO  ,FRIC     ,OLDWSFC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        , GM1      ,ZPI       ,ZPI4GM2,EPSMIN
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOWSHAL  , ONLY : TFAK     ,INDEP    ,CINV
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peakfri.intfb.h"
#include "meansqs_lf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: HALP

      INTEGER(KIND=JWIM) :: IJ, K, M

      ! FREQUENCY INDEXES FOR THE PHILLIPS RANGE 
      INTEGER(KIND=JWIM), PARAMETER :: IPHS = 1 + INT(LOG(1.3_JWRB)/LOG(FRATIO))
      INTEGER(KIND=JWIM), PARAMETER :: IPHE = 1 + INT(LOG(3.0_JWRB)/LOG(FRATIO))

      INTEGER(KIND=JWIM) :: MS, MM, ME
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MMAX

      REAL(KIND=JWRB) :: RF, WFR
      REAL(KIND=JWRB) :: COEF, WS, CHECKTA
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE+IPHE) :: DFRE, F5DFRE 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: F1DMAX
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XMSS 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: F1DWS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FLWS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HALPHAP',0,ZHOOK_HANDLE)

      RF = 0.5_JWRB*(FRATIO - 1.0_JWRB/FRATIO)
      DO M = 1, NFRE
         DFRE(M) = FR(M)*RF
         F5DFRE(M) = FR5(M)*DFRE(M)
      ENDDO
      DO M = NFRE+1, NFRE+IPHE
         DFRE(M) = FR(NFRE)*FRATIO**(M-NFRE)*RF
         F5DFRE(M) = FR5(NFRE)*DFRE(M)
      ENDDO

!     Find windsea spectrum
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

      MMAX(:) = NFRE
      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            WS = 0.5_JWRB + SIGN(0.5_JWRB, (CHECKTA-1.0_JWRB) )
            FLWS(IJ,K,M) = WS*FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      ! Find peak of windsea 1d spectrum
      CALL PEAKFRI (FLWS, IJS, IJL, MMAX, F1DMAX, F1DWS)

      ! Find the Phillips parameter by weighting averaging its value over the Phillips range (see above)
      ! if it is inside the discretised frequency range, otherwise use its relation to the mean square slope
      DO IJ = IJS, IJL
        ALPHAP(IJ) = 0.0_JWRB
        MS = MIN(MMAX(IJ) + IPHS, NFRE)
        MM = MIN(MMAX(IJ) + IPHE, NFRE) 
        ME = MMAX(IJ) + IPHE 

        IF ( ME <= NFRE ) THEN
          WFR = 0.0_JWRB
          DO M = MS, MM 
            WFR = WFR + DFRE(M)
            ALPHAP(IJ) = ALPHAP(IJ) + F5DFRE(M)*F1DWS(IJ,M)
          ENDDO
          ALPHAP(IJ) = ZPI4GM2*ALPHAP(IJ) / WFR
        ELSE IF (MMAX(IJ) < NFRE ) THEN
          CALL MEANSQS_LF(NFRE, IJ, IJ, FL1(IJ:IJ,:,:), XMSS(IJ))
          ALPHAP(IJ) = XMSS(IJ) /LOG(FR(NFRE)/FR(MMAX(IJ)))
        ELSE
          ALPHAP(IJ) = ZPI4GM2*FR5(NFRE)*F1DWS(IJ,NFRE)
        ENDIF
      ENDDO

!     1/2 ALPHAP:
      HALP(:) = 0.5_JWRB*MIN(ALPHAP(:), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
