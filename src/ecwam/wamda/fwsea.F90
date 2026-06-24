      SUBROUTINE FWSEA (KIJS, KIJL, FL1, DEPTH,   &
     &                  EMN, ETOI, USMO, THMO ,   &
     &                  EWFG, FMWFG)

! ----------------------------------------------------------------------

!**** *FWSEA* - ANALYSE THE MODEL SPECTRA TO PROVIDE THE WIND-SEA PART
!               OF THE SPECTRUM, ITS ENERGY, ITS MEAN FREQUENCY.

!      P.LIONELLO     ECMWF    FEBRUARY 1989
!        (MODIFICATION OF A CODE BY P.JANSSEN
!         AND P.LIONELLO - SUMMER '87 - )
!      J.BIDLOT       ECMWF    JANUARY 1997 : CORRECTION FOR INITIAL
!                                             WAVE DIRECTION NOT BEING 0

!     J.BIDLOT       ECMWF    MARCH 2001 : SIMPLIFY THE SEACH BY USING
!                                          THE SAME CRITERIA AS IN
!                                          SEPWISW

!      METHOD.
!      -------

!      A SPECTRAL COMPONENT IS WINDSEA IF
!      1.2 * 28 * U*/C COS(THETA-THMO) > 1

!      WIND SEA ENERGY AND MEAN FREQUENCY ARE COMPUTED BY THE EXTERNAL
!      WSMFEN

!     *KIJS*    FIRST INDEX IN BLOCK.
!     *KIJL*    LAST  INDEX IN BLOCK.
!     *FL1*     MODEL SPECTRUM (FIRST GUESS)
!     *DEPTH*   WATER DEPTH
!     *ETOI*    TOTAL ENERGY FROM O.I. (ANALYSIS)
!     *USMO*    USTAR (FIRST GUESS)
!     *THMO*    WIND DIRECTION
!     *EWFG*    FIRST GUESS WIND-SEA ENERGY
!     *FMWFG*   WINDSEA MEAN FREQUENCY (FIRST GUESS)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH       ,C      ,FRIC     ,OLDWSFC
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wsmfen.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMN, ETOI, USMO, THMO
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EWFG, FMWFG


      INTEGER(KIND=JWIM) :: IJ, K, M, NWP

      REAL(KIND=JWRB) :: THRESHOLD, CHECKTA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CM
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: FSEA

      LOGICAL :: LLFOUND

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FWSEA',0,ZHOOK_HANDLE)

!*   1. PARAMETERS AND INITIALIZATION.
!       ------------------------------

      NWP=0 
      DO IJ=KIJS,KIJL 
        EWFG(IJ)  = -999.0_JWRB
        FMWFG(IJ) = -9999.0_JWRB
      ENDDO


!*    2. LOOP ON THE POINTS INSIDE THE BLOCK.
!        ------------------------------------

      DO M=1,NFRE
        CM(M) = OLDWSFC*FRIC/C(M)
      ENDDO

      DO IJ=KIJS,KIJL
        LLFOUND=.FALSE.
!       SKIP THE LAND POINTS.
        IF (EMN(IJ) > 0.01_JWRB .AND. ETOI(IJ) > 0.01_JWRB) THEN
          DO K=1,NANG
            THRESHOLD = USMO(IJ)*COS(TH(K)-THMO(IJ))
            DO M=1,NFRE
              CHECKTA = CM(M)*THRESHOLD
              IF (CHECKTA < 1.0_JWRB) THEN
                FSEA(K,M) = 0.0_JWRB
              ELSE
                FSEA(K,M) = FL1(IJ,K,M)
                LLFOUND = .TRUE.
              ENDIF
            ENDDO
          ENDDO

!         COMPUTATION OF THE WINDSEA ENERGY WINDSEA MEAN FREQUENCY.
          IF(LLFOUND) THEN
            CALL WSMFEN (FSEA,EWFG(IJ),FMWFG(IJ),USMO(IJ),DEPTH(IJ))
            NWP = NWP + 1
          ENDIF
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FWSEA',1,ZHOOK_HANDLE)

      END SUBROUTINE FWSEA
