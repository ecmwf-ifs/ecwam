      SUBROUTINE SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   NOVEMBER 2017

!*    PURPOSE.
!     --------
!     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
!     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

!**   INTERFACE.
!     ----------
!     *CALL* *SDEPTHLIM((KIJS, KIJL, EMAXDPT, FL1)
!          *KIJS*   - LOCAL INDEX OF FIRST GRIDPOINT
!          *KIJL*   - LOCAL  INDEX OF LAST GRIDPOIN
!          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
!          *FL1*    - SPECTRUM.


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "semean.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EM
      LOGICAL :: LLEPSMIN

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDEPTHLIM',0,ZHOOK_HANDLE)

      LLEPSMIN=.TRUE.
      CALL SEMEAN (FL1, KIJS, KIJL, EM, LLEPSMIN)

      DO IJ=KIJS,KIJL
        EM(IJ)=MIN(EMAXDPT(IJ)/EM(IJ), 1.0_JWRB)
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FL1(IJ,K,M) = MAX(FL1(IJ,K,M)*EM(IJ),EPSMIN) 
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDEPTHLIM',1,ZHOOK_HANDLE)

      END SUBROUTINE SDEPTHLIM
