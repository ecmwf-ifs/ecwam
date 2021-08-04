SUBROUTINE HALPHAP(KIJS, KIJL, USTAR, UDIR, FL1, HALP)

! ----------------------------------------------------------------------

!**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER in the wind direction only


!**   INTERFACE.
!     ----------

!       *CALL* *HALPHAP(KIJS, KIJL, UDIR, FL1, HALP)
!          *KIJS* - INDEX OF FIRST GRIDPOINT
!          *KIJL* - INDEX OF LAST GRIDPOINT
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
      USE YOWSHAL  , ONLY : CINV
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peakfri.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: HALP

      INTEGER(KIND=JWIM) :: IJ, K, M

      ! FREQUENCY INDEXES FOR THE PHILLIPS RANGE 
      INTEGER(KIND=JWIM), PARAMETER :: IPHS = 1 + INT(LOG(1.3_JWRB)/LOG(FRATIO))
      INTEGER(KIND=JWIM), PARAMETER :: IPHE = 1 + INT(LOG(3.0_JWRB)/LOG(FRATIO))

      INTEGER(KIND=JWIM) :: MS, MM, ME
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: MMAX

      REAL(KIND=JWRB) :: RF, WFR
      REAL(KIND=JWRB) :: COEF, WS, CHECKTA
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE+IPHE) :: DFRE, F5DFRE 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1DMAX
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPHAP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: F1DWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLWS

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
        DO IJ = KIJS, KIJL
          XINVWVAGE(IJ,M)=USTAR(IJ)*CINV(IJ,M)
        ENDDO
      ENDDO

      DO K = 1, NANG
        DO IJ = KIJS, KIJL
          DIRCOEF(IJ,K)=COEF*COS(TH(K)-UDIR(IJ))
        ENDDO
      ENDDO

      MMAX(:) = NFRE
      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = KIJS, KIJL
            CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            WS = 0.5_JWRB + SIGN(0.5_JWRB, (CHECKTA-1.0_JWRB) )
            FLWS(IJ,K,M) = WS*FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      ! Find peak of windsea 1d spectrum
      CALL PEAKFRI (KIJS, KIJL, FLWS, MMAX, F1DMAX, F1DWS)

      ! Find the Phillips parameter by weighting averaging its value over the Phillips range (see above)
      DO IJ = KIJS, KIJL
        ALPHAP(IJ) = 0.0_JWRB
        MS = MIN(MMAX(IJ) + IPHS, NFRE)
        MM = MIN(MMAX(IJ) + IPHE, NFRE) 
        ME = MMAX(IJ) + IPHE 
        WFR = 0.0_JWRB
        DO M = MS, MM 
          WFR = WFR + DFRE(M)
          ALPHAP(IJ) = ALPHAP(IJ) + F5DFRE(M)*F1DWS(IJ,M)
        ENDDO
        ! extension above FR(NFRE) with f**-5 tail
        DO M = NFRE+1, ME
          WFR = WFR + DFRE(M)
          ALPHAP(IJ) = ALPHAP(IJ) + F5DFRE(M)*F1DWS(IJ,NFRE)
        ENDDO
        ALPHAP(IJ) = ZPI4GM2*ALPHAP(IJ) / WFR
      ENDDO

!     1/2 ALPHAP:
      HALP(:) = 0.5_JWRB*MIN(ALPHAP(:), ALPHAPMAX)

IF (LHOOK) CALL DR_HOOK('HALPHAP',1,ZHOOK_HANDLE)

END SUBROUTINE HALPHAP
