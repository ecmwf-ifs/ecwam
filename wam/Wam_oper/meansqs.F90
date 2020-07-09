      SUBROUTINE MEANSQS(IJS, IJL, F, USTAR, THW, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE.

!     P.A.E.M. JANSSEN
!     J. BIDLOT  ECMWF  FEBRUARY 1996  MESSAGE PASSING

!     July 2020, add gravity capillary contribution


!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (IJS, IJL, F, USTAR, THW, XMSS)*
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *F*   - SPECTRUM.
!              *USTAR* - NEW FRICTION VELOCITY IN M/S (INPUT).
!              *THW*  - WIND DIRECTION
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

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,DFIM_SIM ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "halphap.intfb.h"
#include "meansqs_gc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THW 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XMSS 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JWRB), PARAMETER :: XLAMBDAC=0.0628_JWRB

      REAL(KIND=JWRB) :: FS, XKC, FC, CONST1, CONST2, CP, XI, ALPHAP, CONST3 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE_ODD) :: FD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) ::  TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: HALP, FRGC

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MEANSQS',0,ZHOOK_HANDLE)

!*    1. COMPUTE THE GRAVITY-CAPILLARY CONTRIBUTION TO MSS 
!        -------------------------------------------------

!     COMPUTE THE PHILLIPS PARAMETER
      CALL HALPHAP(IJS, IJL, F, THW, HALP)

!     GRAVITY-CAPILLARY CONTRIBUTION TO MSS
      CALL MEANSQS_GC(IJS, IJL, HALP, USTAR, XMSS, FRGC)


!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE_ODD
          FD(M) = DFIM_SIM(M)*(ZPIFR(M))**4/G**2
        ENDDO

        DO M=1,NFRE_ODD
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

        DO M=1,NFRE_ODD
          DO IJ=IJS,IJL
            TEMP1(IJ) = DFIM_SIM(M)*TFAK(INDEP(IJ),M)**2
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

!*    3. ADD TAIL CORRECTION TO MEAN SQUARE SLOPE.
!        ------------------------------------------

!! not applied !!
!!      FS     = FR(NFRE_ODD)      
!!      XKC    = ZPI/XLAMBDAC
!!      FC     = SQRT(G*XKC+SURFT*XKC**3)/ZPI
!!      CONST2 = ZPI**4*FS**5/G**2
!!      DO IJ=IJS,IJL      
!!        CONST1    = 0.0_JWRB*LOG(FC/FS)  
!!        CP        = G/(ZPI*FM(IJ))
!!        XI        = CP/USNEW(IJ)
!!        ALPHAP    = MAX(0.21_JWRB/XI,0.0040_JWRB)
!!        CONST3    = CONST2*DELTH*TEMP2(IJ)
!!        XMSS(IJ) = XMSS(IJ)+CONST1*CONST3
!!      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS
