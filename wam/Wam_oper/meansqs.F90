      SUBROUTINE MEANSQS(XKMSS, IJS, IJL, F, U10, USTAR, UDIR, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE UP TO WAVE NUMBER XKMSS.

!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (XKMSS, IJS, IJL, F, U10, USTAR, UDIR, XMSS)*
!              *XKMSS* - WAVE NUMBER CUT OFF
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *F*   - SPECTRUM.
!              *U10* - 10m wind speed
!              *USTAR* - NEW FRICTION VELOCITY IN M/S (INPUT).
!              *UDIR* - WIND SPEED DIRECTION
!              *XMSS* - MEAN SQUARE SLOPE (OUTPUT).

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI   , ZPI4GM2
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,DFIM  , DELTH, FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "halphap.intfb.h"
#include "meansqs_gc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XMSS 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      INTEGER(KIND=JWIM) :: IJ, M, K, NFRE_MSS, NFRE_EFF

      REAL(KIND=JWRB) :: XLOGFS, FCUT
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XMSS_TAIL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: HALP, FRGC

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MEANSQS',0,ZHOOK_HANDLE)

!*    1. COMPUTE THE GRAVITY-CAPILLARY CONTRIBUTION TO MSS 
!        -------------------------------------------------

!     COMPUTE THE PHILLIPS PARAMETER
      CALL HALPHAP(IJS, IJL, UDIR, F, HALP)

!     GRAVITY-CAPILLARY CONTRIBUTION TO MSS
      CALL MEANSQS_GC(XKMSS, IJS, IJL, HALP, U10, USTAR, XMSS, FRGC)

      FCUT = SQRT(G*XKMSS)/ZPI
      NFRE_MSS = INT(LOG(FCUT/FR(1))/LOG(FRATIO))+1 
      NFRE_EFF = MIN(NFRE,NFRE_MSS)

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE_EFF
          FD(M) = DFIM(M)*(ZPIFR(M))**4*GM1**2
        ENDDO

        DO M=1,NFRE_EFF
          DO IJ=IJS,IJL
            TEMP2(IJ) = 0.0_JWRB
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            XMSS(IJ) = XMSS(IJ)+FD(M)*TEMP2(IJ)
          ENDDO
        ENDDO
!SHALLOW
      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE_EFF
          DO IJ=IJS,IJL
            TEMP1(IJ) = DFIM(M)*TFAK(INDEP(IJ),M)**2
            TEMP2(IJ) = 0.0_JWRB
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            XMSS(IJ) = XMSS(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO
      ENDIF
!SHALLOW

!*    3. ADD TAIL CORRECTION TO MEAN SQUARE SLOPE (between FR(NFRE_EFF) and FRGC).
!        ------------------------------------------

      XLOGFS = LOG(FR(NFRE_EFF))
      DO IJ=IJS,IJL
        XMSS_TAIL(IJ) = 2.0_JWRB*HALP(IJ)*MAX((LOG(MIN(FRGC(IJ),FCUT)) - XLOGFS),0.0_JWRB)
        XMSS(IJ) = XMSS(IJ) + XMSS_TAIL(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS
