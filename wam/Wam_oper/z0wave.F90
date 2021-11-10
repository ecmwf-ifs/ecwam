      SUBROUTINE Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B)

! ----------------------------------------------------------------------

!**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.

!*    PURPOSE.
!     --------

!       COMPUTE ROUGHNESS LENGTH. 

!**   INTERFACE.
!     ----------

!       *CALL* *Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0)
!          *KIJS* - INDEX OF FIRST GRIDPOINT.
!          *KIJL* - INDEX OF LAST GRIDPOINT.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *UTOP* - WIND SPEED.
!          *Z0*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *Z0B*  - BACKGROUND ROUGHNESS LENGTH.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK
      USE YOWPCONS , ONLY : GM1
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWTABL  , ONLY : EPS1

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "chnkmin.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL),INTENT(IN)  ::  US, TAUW, UTOP
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL),INTENT(OUT) ::  Z0, Z0B


      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: UST2, UST3, ARG
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPHAOG

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('Z0WAVE',0,ZHOOK_HANDLE)

      IF (LLCAPCHNK) THEN
        DO IJ=KIJS,KIJL
          ALPHAOG(IJ)= CHNKMIN(UTOP(IJ))*GM1
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          ALPHAOG(IJ)= ALPHA*GM1
        ENDDO
      ENDIF

      DO IJ=KIJS,KIJL
        UST2 = US(IJ)**2
        UST3 = US(IJ)**3
        ARG = MAX(UST2-TAUW(IJ),EPS1)
        Z0(IJ) = ALPHAOG(IJ)*UST3/SQRT(ARG)
        Z0B(IJ) = ALPHAOG(IJ)*UST2
      ENDDO

      IF (LHOOK) CALL DR_HOOK('Z0WAVE',1,ZHOOK_HANDLE)

      END SUBROUTINE Z0WAVE
